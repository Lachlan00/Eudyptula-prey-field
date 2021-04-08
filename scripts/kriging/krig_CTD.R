#' 3D-CTD Kriging
#' 
#' Authors: Sven Gastauer & Lachlan Phillips
#' ACE CRC, 2018 (Hobart, Tasmania, Australia)
#' MQ Uni, 2020 (Syndey, New South Wales, Australia)
#' 
#' Modified by Lachlan Phillips 2019-Dec
#' Macquarie University - Sydney, NSW
#' Original code can be found in this directory
#' "sven/R/CTD3DKriging.R"
#'
#' Load CTD data, create a grid and data based variograms to krige the data in 3 Dimensions.
#' Needed inputs are the rdepth reolution at which the data should enter the Geostats database,
#' the resolution of the projected kriging map and the CTD variable which should be kriged
#' 
#' 1 Preparation
#'  1.1 Libraries
#'  1.2 Directories
#'  1.3 Helper Functions - deprecated
#'  1.4 Model parameters
#'
#' 2 Load the data
#'  2.1 Log data - deprecated
#'  2.2 CTD Data
#'  
#' 3 Kriging and plotting   
#'  3.1 Geostats
#'  3.2 Plotting
#'

#################
# 1 Preparation #
#################

#################
# 1.1 Libraries #
#################
library(rgl)
library(plot3D)
library(plot3Drgl)
library(ggplot2)
library(reshape2)
library(viridis)
library(cmocean)
library(beepr)
library(RGeostats)
library(RJSONIO)

###################
# 1.2 Directories #
###################
setwd("/Users/lachlanphillips/Development/PhD/repos/Eudyptula")

########################
# 1.3 Helper functions #
########################
# Lachlan's CTD functions
source('scripts/CTD_functions.R')
source('scripts/utilities.R')

###################################
# 1.4 Define the model parameters #
###################################
# Define which variable should be kriged
vars <- c("temperature","salinity") # 'all' = ["temperature","sound_velocity","salinity","density"]
# Define dr (depth resolution) the resolution of the depth variable used as
# an input into the model. Note this is not the final grid resolution
dr=1
# Define the grid size in 3D 
# [Lon resolution in dec deg, Lat resolution in dec deg, Depth resolution in m]
gridsize=c(0.2,0.2,1) #c(.002,.002,1) # .002 is good, .1 for NM
polygon_fn <- 'kriging/assets/kriging_field/survey-coast-straight_poly.RDS'
# neighbourhood paramters
xy_res=2 # rounding factor for lat/lon aggregation
nmini=3 # Minimum number of points in the neighborhood 
nmaxi=10 # Maximum number of points in the neighborhood
radius=c(50,50,20) # Search radius (in km for UTM and meters for depth)
flag.sector=TRUE # use sector searches?
nsect=3 # direction sectors to use
nsmax=4 # maximum number of points per sector
flag.continuous=FALSE # continous moving neighbourhood?
flag.linked=FALSE # link variabkle drift?
# output directories
model_out_dir <- "./kriging/output/CTD/models/"
plots_out_dir <- "./kriging/output/CTD/plots/"
paramter_json_fn <- "./kriging/output/CTD/model_settings.json"
# set plotting colours
make_plots <- FALSE # better plots made in alternative function
colour_maps <- list(temperature = cmocean('thermal')(100),
                    salinity = cmocean('haline')(100),
                    soundVelocity = cmocean('speed')(100),
                    density = cmocean('dense')(100))

# Debug mode (saves into debug folder with parametrs saved)
DEBUG <- FALSE
debug_var = 'temperature'
debug_survey = '2019_S2'

##########################
# Log parameter settings #
##########################
if (!DEBUG){
  params <- list(vars=vars,
                 dr=dr,
                 gridsize=gridsize,
                 polygon_fn=polygon_fn,
                 xy_res=xy_res,
                 nmini=nmini,nmaxi=nmaxi,radius=radius,
                 flag.sector=flag.sector,nsect=nsect,nsmax=nsmax,
                 flag.continuous=flag.continuous,flag.linked=flag.linked)
  write(toJSON(params), paramter_json_fn)
}

##################
# Debug settings #
##################
# make debug directories
if (DEBUG){
  if (flag.sector){
    sector_str <- paste0('__+flag.sector_nsect-',nsect,'_nsmax-',nsmax)
  } else {
    sector_str <- ''
  }
  if (flag.continuous){
    cont_str <- '__+flag.continuous'
  } else {
    cont_str <- ''
  }
  if (flag.linked){
    link_str <- '__+flag.linked'
  } else {
    link_str <- ''
  }
  debug_dir <- paste0("./kriging/output/debug/param_results/CTD/nmini-",
                      nmini,"__nmaxi-",nmaxi,"__radius-",paste(radius,collapse='-'),"__xyres-",
                      xy_res,sector_str,cont_str,link_str)
  model_out_dir <- paste0(debug_dir,'/models/')
  plots_out_dir <- paste0(debug_dir,'/plots/')
  if (!dir.exists(debug_dir)){
    dir.create(debug_dir)
    dir.create(model_out_dir)
    dir.create(plots_out_dir)
  } else {
    message('Warning: Settings already run...')
    beep(sound='fanfare')
  }
}

################
# 2.2 CTD Data #
################
# **original code LARGELY REPLACED BY `cast.reader()`** #
CTD_data <- cast.reader(drop_non_transect=TRUE)
# filter to the survey set
CTD_data <- cast.survey.filter(CTD_data)
# split out the meta and data frame
CTD_df <- CTD_data[[1]]
# expand transect data
CTD_df <- expand.transect.ID(CTD_df)
# add in lat lon from CTD meta data
CTD_df <- CTD.add.grid.coords(CTD_df)
# Round the depth to the input parameter dr
CTD_df$rdepth <- round(CTD_df$depth/dr,0)*dr

# aggregate data
CTD_df <- aggregate(CTD_df, by=list(CTD_df$rdepth, 
                                    paste0(CTD_df$survey_id,'_',CTD_df$transect_id),
                                    round(CTD_df$lat, xy_res), round(CTD_df$lon, xy_res)),
                    FUN="mean", na.action=na.pass)
# sub back in survey id
CTD_df$survey_id <- substr(CTD_df$Group.2, 1, 7)
# chnage sound velocity variable name so it doesn't mess with file naming system
colnames(CTD_df)[colnames(CTD_df) == "sound_velocity"] <- 'soundVelocity'

#######################
# 3.1 Initialise grid #
#######################
# read polygon
poly <- readRDS(polygon_fn)
# Correct polygon Lon and Lat to be in UTM (km)
poly <- lonlat2UTM.poly(poly)
poly@sets[[1]]$x = poly@sets[[1]]$x/1000
poly@sets[[1]]$y = poly@sets[[1]]$y/1000
# setup intiial grid so consistent between each run
ccc <- subset(CTD_df[,c("lat","lon","rdepth","temperature")])
# get lon and lat in UTM (km)
ccc <- lonlat2UTM.df(ccc)
ccc$lat =  ccc$lat/1000
ccc$lon =  ccc$lon/1000
# add one more point to extend grid to edge of polygon
ccc <- rbind(ccc, data.frame(lat=-4032, lon=min(poly@sets[[1]]$x), rdepth=0, temperature=16))
# create database and define dimensions
db.ctd <- db.create(ccc)
db.ctd <- db.locate(db.ctd, c("lon","lat","rdepth"),"x")
db.ctd <- db.locate(db.ctd, c("temperature"),"z")
#initialise a grid
db.grid <- db.grid.init(db.ctd, dcell=gridsize)
# make grid from polygon
db.grid <- db.polygon(db.grid, poly)

############################
# 3.2 Kriging and plotting #
############################
if (vars=='all'){
  vars <- c("temperature","soundVelocity","salinity","density")
}
# reduce
if (DEBUG){
  vars = debug_var
  CTD_df <- CTD_df[CTD_df$survey_id == debug_survey,]
  message(paste('DEBUG mode TRUE - kriging',debug_var,'-',debug_survey))
}
for(var in vars){
  message('\n',Sys.time(),": Processing ",var)
  # loop through the different surveys
  for(survey in unique(CTD_df$survey_id)){
    message('\t',Sys.time(),": Processing ",survey)
    # select survey subset
    ctd_sub <- CTD_df[CTD_df$survey_id==survey,]
    # create new dataframe only containing relevant information
    ccc <- subset(ctd_sub[,c("lat","lon","rdepth",var)])
    # convert to UTM m and then to km (better for kriging)
    ccc <- lonlat2UTM.df(ccc)
    ccc$lat =  ccc$lat/1000
    ccc$lon =  ccc$lon/1000
    ##############################################
    # 3.1 Geostats
    ##############################################
    # create database and define dimensions
    db.ctd <- db.create(ccc)
    db.ctd <- db.locate(db.ctd, c("lon","lat","rdepth"),"x")
    db.ctd <- db.locate(db.ctd, c(var),"z")
    
    # Calculate experimental variogram
    vari <- vario.calc(db.ctd)
    plot(vari)
    #Fit variogram mode
    vg.mod <- model.auto(vari)
    # Define neighbourhood
    neimov <- neigh.create(ndim=3,type=2,nmini=nmini,nmaxi=nmaxi,radius=radius,
                           flag.sector=flag.sector,nsect=nsect,nsmax=nsmax,
                           flag.continuous=flag.continuous)
    # Run the kriging
    message('\t\t',Sys.time(),': Starting kriging...')
    kres <- kriging(dbin=db.ctd, dbout=db.grid, model=vg.mod, 
                    neigh=neimov, flag.linked=flag.linked)
    message('\t\t',Sys.time(),': Saving kriging...')
    # convert back to UTM (m)
    kres <- kres_UTMkm2UTMm(kres)
    # save database
    saveRDS(kres, file=paste0(model_out_dir,"3Dkrig_CTDdb_",var,"_",survey,".rds"))
    # Convert kriging results into dataframe
    t.df <- as.data.frame(kres@items)#[kres@items$Polygon==TRUE,])
    names(t.df)<- c('rank', 'lon','lat','depth','polygon',var,'stdev')
    # Convert Lon lat into degrees
    t.df <- UTM2lonlat.df(t.df)

    # add survey information
    t.df$survey <- survey
    #Set all data outside of Polygon to be NA
    t.df[,var][t.df$polygon==FALSE] <- NA
    # dataframes not used so probably do not need anymore
    # saveRDS(t.df, file=paste0(model_out_dir,"3Dkrig_CTDdf_",var,"_",survey,".rds"))

    # # Make high contrast cast stations
    # t.df <- readRDS(paste0("./kriging/output/CTD/models/3Dkrig_CTDdf_",var,"_",survey,".RDS"))
    # cast_idx <- krig.get.cast.index(t.df)
    # t.df[cast_idx,var] <- 1e4

    #####################
    # 3.2 Visualisation #
    #####################
    message('\t\t',Sys.time(),': Plotting...')
    if (make_plots){
      # temp fix to fix vector memory errors
      t.df <- lonlat2UTM.df(t.df)
      #Set up dimension
      xrange <- unique(t.df$lon)
      yrange <- unique(t.df$lat)
      zrange <- -unique(t.df$depth)
      x <- xrange; y <- yrange; z <- zrange
      #Create Mesh
      M <- mesh(xrange, yrange, zrange)
      #Create 3D data array
      p <- array(t.df[,var], dim=c(length(xrange),length(yrange), length(zrange)))

      # Make open 3D plot
      png(paste0(plots_out_dir,"3dKrig_",var,"_",survey,"_open.png"),
          width=1000, height=850)
      # add margin
      # par(mar=c(5.1, 7.1, 4.1, 2.1))

      slice3D(x, y, z, colvar = p,
              col = colour_maps[[var]],
              NAcol="transparent",
              xs=NULL,
              ys=yrange[c(30,50,70,90)],
              zs=c(zrange[c(1,20,40,60,80,100)]),
              d = 1,
              aspect=c(1,1,0.1),
              xlab="Longitude",
              ylab="Latitude",
              zlab="Depth",
              clab=var,
              ticktype="detailed")

      dev.off()

      # Make closed 3D plot
      png(paste0(plots_out_dir,"3dKrig_",var,"_",survey,"_close.png"),
          width=1000, height=850)

      slice3D(x, y, z, colvar = p,
              col = colour_maps[[var]],
              NAcol="transparent",
              xs=NULL,
              ys = yrange[c(1,30,60,90)],
              zs=c(zrange[c(1,20,40,60,80,100)]),
              d = 1,
              xlab="Longitude",
              ylab="Latitude",
              zlab="Depth",
              clab=var,
              ticktype="detailed")

      dev.off()

      # cross section image
      t <- t.df[t.df$polygon==TRUE,]

      xsection <- ggplot(data=t[t$depth %in% c(6,10,20,40),],
             aes(x=lon, y=lat,fill=t[t$depth %in% c(6,10,20,40),var])) +
        geom_raster() +
        scale_fill_gradientn(colours=colour_maps[[var]]) +
        labs(fill=col_lab) +
        facet_wrap(~depth)

      ggsave(paste0(plots_out_dir,"3Dkrig_xsection_",var,"_",survey,".png"),
             xsection, width=5, height=5, dpi=200)
      # end processing
    }
  }
  message(Sys.time(),": Completed ",var)
}

# Alert script end
beep()
