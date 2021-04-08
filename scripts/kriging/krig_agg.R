#' 3D-prey-field Kriging
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
#' Load acoustic data, create a grid and data based variograms to krige the data in 3 Dimensions.
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
library(RGeostats)
library(ggplot2)
library(ggpubr)
library(ggmap)
library(reshape2)
library(plot3Drgl)
library(rgl)
library(plot3D)
library(cmocean)
library(beepr)
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
# School aggregation data directory
agg_dir <- '/Volumes/LP_MstrData/master-data/survey/acoustics/acoustics_processed/'
# Define which variable should be kriged
vars <- 'Sv_mean' # all = c('Sv_mean','NASC')
# Define dr (depth resolution) the resolution of the depth variable used as
# an input into the model. Note this is not the final grid resolution
dr=1
# Define the grid size in 3D 
gridsize=c(0.2,0.2,1) # c(km, km, m)
polygon_fn <- 'kriging/assets/kriging_field/survey-coast-straight_poly.rds'
# neighbourhood paramters
xy_res=2 # rounding factor for lat/lon aggregation
nmini=3 # Minimum number of points in the neighborhood 
nmaxi=10 # Maximum number of points in the neighborhood
radius=c(50,50,20) # Search radius (in km for UTM and meters for depth)
flag.sector=FALSE # use sector searches?
nsect=2 # direction sectors to use
nsmax=12 # maximum number of points per sector
flag.continuous=TRUE # continous moving neighbourhood?
flag.linked=FALSE # link variable drift?
linear.transform=TRUE # Transform to linear space?
normalise.data=TRUE # Transform data to be normal prior to kriging
# output directories
model_out_dir <- "./kriging/output/agg/models/"
plots_out_dir <- "./kriging/output/agg/plots/"
paramter_json_fn <- "./kriging/output/agg/model_settings.json"
# set plotting colours
make_plots <- TRUE # better plots made in alternative function
cmap <- rev(cmocean('deep')(100))
# Checks
run.checks=FALSE # Run model checks  (not finished)
 
# Debug mode (saves into debug folder with parametrs saved)
DEBUG <- FALSE
debug_var = 'Sv_mean'
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
                 flag.continuous=flag.continuous,flag.linked=flag.linked,
                 linear.transform=linear.transform,
                 normalise.data=normalise.data)
  write(toJSON(params), paramter_json_fn)
}

##################
# Debug settings #
##################
# make debug directories
if (DEBUG){
  sector_str  <- ifelse(flag.sector, paste0('__+flag.sector_nsect-',nsect,'_nsmax-',nsmax), '')
  cont_str <- ifelse(flag.continuous, '__+flag.continuous', '')
  link_str <- ifelse(flag.linked, '__+flag.linked', '')
  linear_str <- ifelse(linear.transform, '__+linearTrans', '')
  norm_str <- ifelse(normalise.data, '__+normalised', '')
  debug_dir <- paste0("./kriging/output/debug/param_results/agg/nmini-",
                      nmini,"__nmaxi-",nmaxi,"__radius-",paste(radius,collapse='-'),"__xyres-",
                      xy_res,sector_str,cont_str,link_str,linear_str,norm_str)
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

#################################
# 2.1 CTD Data (for grid setup) #
#################################
# CTD data is called so kriging grid is 
# consistent between CTD and preyfield 
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
# change sound velocity variable name so it doesn't mess with file naming system
colnames(CTD_df)[colnames(CTD_df) == "sound_velocity"] <- 'soundVelocity'

#######################
# 2.2 Initialise grid #
#######################
# Initialise a consistent grid across all surveys
# Use same grid as CTD data
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
# initialise a grid
db.grid <- db.grid.init(db.ctd, dcell=gridsize)
# make grid from polygon
db.grid <- db.polygon(db.grid, poly)

#######################################
# 3.1 Load in school aggregation data #
#######################################
if (vars=='all'){
  vars <- c('Sv_mean', 'NASC')
}
# list all files
file.ls <- list.files(agg_dir, recursive = T)
# get aggregations and gps
# doing some grep to filter out secondary acoustic sets that need to be untangaled
agg.ls <- file.ls[grepl('*[S][1-2]_schools.csv', file.ls)]
gps.ls <- file.ls[grepl('*[S][1-2].gps.csv', gsub('_', '', file.ls))]
# because some datasets don't have separate gps files we need to get survey ids for each
agg.ls.surveys <- substr(agg.ls, 1, 7)
gps.ls.surveys <- substr(gps.ls, 1, 7)

# NOTE!!
# LATER - make sure acoustic transections and CTD krig have proper overlap
# When this is done make sure to use secondary acoustic sets filtered above

# Now load in the data
agg.ls <- lapply(paste0(agg_dir, agg.ls), read.csv)
# for each data frame attach the survey id and then filter to the data
# we want to krig
keeps <- c('survey_id','dtUTC','Depth_mean','Lon_S','Lat_S',vars)
for (i in 1:length(agg.ls)){
  # add survey information
  agg.ls[[i]]$survey_id <- agg.ls.surveys[i]
  # calc datetimes
  agg.ls[[i]]$dtUTC <- as.POSIXct(paste(agg.ls[[i]]$Date_S, agg.ls[[i]]$Time_S), 
                                  format='%Y%m%d %H:%M:%OS', tz='UTC')
  # fiter bad GPS data
  agg.ls[[i]] <- agg.ls[[i]][agg.ls[[i]]$Lon_S <= 360,]
  agg.ls[[i]] <- agg.ls[[i]][agg.ls[[i]]$Lat_S <= 90,]
  # filter bad values
  agg.ls[[i]] <- agg.ls[[i]][agg.ls[[i]]$NASC < 50000,]
  agg.ls[[i]] <- agg.ls[[i]][agg.ls[[i]]$NASC > -0.1,]
  agg.ls[[i]] <- agg.ls[[i]][agg.ls[[i]]$Sv_mean > -500,]
  # restrict data to what we want
  agg.ls[[i]] <- agg.ls[[i]][,keeps]
}
# now merge the data into 1 dataframe
agg.data <- do.call(rbind, agg.ls)
# Round the depth to the input parameter dr
agg.data$Depth_mean <- round(agg.data$Depth_mean/dr,0)*dr
colnames(agg.data)[colnames(agg.data) == 'Depth_mean'] <- 'rdepth'
colnames(agg.data)[colnames(agg.data) == 'Lon_S'] <- 'lon'
colnames(agg.data)[colnames(agg.data) == 'Lat_S'] <- 'lat'
agg.data$dtUTC <- NULL
# aggregate the data
agg.data <- aggregate(agg.data, by=list(agg.data$rdepth, 
                                        agg.data$survey_id,
                                        round(agg.data$lat, xy_res), 
                                        round(agg.data$lon, xy_res)),
                      FUN="mean", na.action=na.pass)
# sub back in survey id
agg.data$survey_id <- agg.data$Group.2
# organise data
agg.data <- agg.data[,5:ncol(agg.data)]

############################
# 4.1 Kriging and plotting #
############################
# reduce vars if DEBUG mode is TRUE
if (DEBUG){
  vars = debug_var
  agg.data <- agg.data[agg.data$survey_id == debug_survey,]
  message(paste('DEBUG mode TRUE - kriging',debug_var,'-',debug_survey))
}
for(var in vars){
  message('\n',Sys.time(),": Processing ",var)
  # loop through the different surveys
  for(survey in unique(agg.data$survey_id)){
    message('\t',Sys.time(),": Processing ",survey)
    # create datavase from data subset
    agg_sub <- agg.data[agg.data$survey_id==survey,]
    # Transform data
    if (linear.transform){
      agg_sub[var] <- (10**(agg_sub[var]/10))
    }
    if (normalise.data){
      normal.trans <- nscore(agg_sub[[var]])
      agg_sub[var] <- normal.trans$nscore
    }
    # create new dataframe only containing relevant information
    ccc <- subset(agg_sub[,c("lat","lon","rdepth",var)])
    # convert to UTM m and then to km (better for kriging)
    ccc <- lonlat2UTM.df(ccc)
    ccc$lat =  ccc$lat/1000
    ccc$lon =  ccc$lon/1000
    ################
    # 4.2 Geostats #
    ################
    # create database and define dimensions
    db.agg <- db.create(ccc)
    db.agg <- db.locate(db.agg, c("lon","lat","rdepth"),"x")
    db.agg <- db.locate(db.agg, c(var),"z")
    # Calculate experimental variogram
    vari <- vario.calc(db.agg)
    plot(vari)
    # Fit variogram mode
    vg.mod <- model.auto(vari,
                         struct = c('Nugget Effect', 'Gaussian', 'Cubic'),
                         flag.noreduce=T)
    # Define neighbourhood
    neimov <- neigh.create(ndim=3,type=2,nmini=nmini,nmaxi=nmaxi,radius=radius,
                           flag.sector=flag.sector,nsect=nsect,nsmax=nsmax,
                           flag.continuous=flag.continuous, flag.aniso=TRUE)
    # # Run checks
    # if (run.checks){
    #   # cross validate neighbourhood
    #   db.data.xv <- xvalid(db=db.agg, model=vg.mod, neigh=neimov, uc=c("1"),
    #                        mean=NA)
    #   # Mean error (should be close to 0) in percent of zone mean
    #   mean(db.extract(db.data.xv, paste0('Xvalid.',var,'.esterr'))) / 963.26
    #   ## 0.0453
    #   # Mean standardized squared error (should be close to 1)
    #   mean(db.extract(db.data.xv, "Xvalid.eggs.stderr")^2)
    #   ## 1.2591
    # }
    # Run the kriging
    message('\t\t',Sys.time(),': Starting kriging...')
    kres <- kriging(dbin=db.agg, dbout=db.grid, model=vg.mod,
                    neigh=neimov, flag.linked=flag.linked)
    message('\t\t',Sys.time(),': Saving kriging...')
    kres <- kres_UTMkm2UTMm(kres)
    # Back transform data
    if (normalise.data){
      # Save transformation information
      normal.trans$transformed.krig.res <- kres@items[paste0('Kriging.',var,'.estim')]
      kres@items[paste0('Kriging.',var,'.estim')] <- backtr.nscore(kres@items[[paste0('Kriging.',var,'.estim')]],
                                                                   normal.trans)
      normal.trans$linear.trans <- FALSE
    }
    # Linear to log space
    if (linear.transform){
      kres@items[paste0('Kriging.',var,'.estim')] <- 10*log10(kres@items[paste0('Kriging.',var,'.estim')])
      normal.trans$linear.trans <- TRUE
    }
    # Save trasnformation information
    if (normalise.data){
      saveRDS(normal.trans, file=paste0(model_out_dir,"transformation_aggdb_",var,"_",survey,".rds"))
    }
    # save database
    saveRDS(kres, file=paste0(model_out_dir,"3Dkrig_aggdb_",var,"_",survey,".rds"))
    # save variogram
    saveRDS(vari, file=paste0(model_out_dir,"vari_aggdb_",var,"_",survey,".rds"))
    # Convert kriging results into dataframe
    t.df <- as.data.frame(kres@items)#[kres@items$Polygon==TRUE,])
    names(t.df)<- c('rank', 'lon','lat','depth','polygon',var,'stdev')
    # Convert Lon lat into degrees
    t.df <- UTM2lonlat.df(t.df)
    # add survey information
    t.df$survey <- survey
    #Set all data outside of Polygon to be NA
    t.df[,var][t.df$polygon==FALSE] <- NA

    #####################
    # 4.3 Visualisation #
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
              col = cmap,
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
              col = cmap,
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

      xsection <- ggplot(data=t[t$depth %in% c(5,10,20,40),],
                         aes(x=lon, y=lat,fill=t[t$depth %in% c(6,10,20,40),var])) +
        geom_raster() +
        scale_fill_gradientn(colours=cmap) +
        labs(fill=var) +
        facet_wrap(~depth)
      ggsave(paste0(plots_out_dir,"3Dkrig_xsection_",var,"_",survey,".png"),
             xsection, width=5, height=5, dpi=200)

    }
    beep()
  }
  message(Sys.time(),": Completed ",var)
}

# Method for 3d sctater plot
# scatterplot3d(agg.data$lon, agg.data$lat, -agg.data$rdepth, angle = 55, 
#               color=map2color(-agg.data$rdepth,rev(cmocean('deep')(200)[50:200])))