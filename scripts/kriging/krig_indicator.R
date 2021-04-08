# Indicator kriging
# This script krigs the probability of acoustic thresholds 
# at different depths for Sv_mean

#################
# 1.1 Libraries #
#################
library(RGeostats)
library(ggplot2)
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
var <- 'Sv_mean' # all = c('Sv_mean','NASC')
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
flag.sector=FALSE # use sector searches?
nsect=4 # direction sectors to use
nsmax=7 # maximum number of points per sector
flag.continuous=TRUE # continous moving neighbourhood?
flag.linked=FALSE # link variabkle drift?
# output directories
model_out_dir <- "./kriging/output/indicator/models/"
plots_out_dir <- "./kriging/output/indicator/plots/"
paramter_json_fn <- "./kriging/output/indicator/model_settings.json"
# set plotting colours
make_plots <- TRUE # better plots made in alternative function
cmap <- 'haline'

# Debug mode (saves into debug folder with parametrs saved)
DEBUG <- FALSE
debug_var = 'Sv_mean'
debug_survey = '2019_S1' # next 2017_S2

##########################
# Log parameter settings #
##########################
if (!DEBUG){
  params <- list(var=var,
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
  debug_dir <- paste0("./kriging/output/debug/param_results/indicator/nmini-",
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
# chnage sound velocity variable name so it doesn't mess with file naming system
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
#initialise a grid
db.grid <- db.grid.init(db.ctd, dcell=gridsize)
# make grid from polygon
db.grid <- db.polygon(db.grid, poly)

#######################################
# 3.1 Load in school aggregation data #
#######################################
# list all files
file.ls <- list.files(agg_dir, recursive = T)
# get aggregations and gps
# doing some grep to filter out secondary acoustic sets that need to be untangaled
agg.ls <- file.ls[grepl('*[S][1-2]_schools.csv', file.ls)]
# gps.ls <- file.ls[grepl('*[S][1-2].gps.csv', gsub('_', '', file.ls))]
# because some datasets don't have separate gps files we need to get survey ids for each
agg.ls.surveys <- substr(agg.ls, 1, 7)
# gps.ls.surveys <- substr(gps.ls, 1, 7)

# NOTE!!
# LATER - make sure acoustic transections and CTD krig have proper overlap
# When this is done make sure to use secondary acoustic sets filtered above

# Now load in the data
agg.ls <- lapply(paste0(agg_dir, agg.ls), read.csv)
# for each data frame attach the survey id and then filter to the data
# we want to krig
keeps <- c('survey_id','dtUTC','Depth_mean','Lon_S','Lat_S','Sv_mean','NASC')
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
# convert coordinates
agg.data <- lonlat2UTM.df(agg.data)
agg.data$lon <- agg.data$lon/1000
agg.data$lat <- agg.data$lat/1000
# reduce to what we want
agg.data_all <- agg.data
keeps <- c('survey_id','rdepth','lon','lat',var)
agg.data <- agg.data[,keeps]

#########################################
# 4.1 Define acoustic threshold cutoffs #
#########################################
# make a databse of all surveys
db.sub <- db.create(agg.data)
db.sub <- db.locate(db.sub,c("lon","lat","rdepth"),"x")
db.sub <- db.locate(db.sub, var, "z")
# set thresholds from quantiles
zcut <- as.numeric(quantile(db.sub[, var])) #c(-80,-60,-50,-40)
my.limits <- limits.create(zcut=zcut[-5], flag.zcut.int = F)
db.sub <- db.indicator(db.sub, my.limits)

############################
# 4.2 Calculate variograms #
############################
# Variograhy
dirvect <- c(0,90)
# Mean annual variogram
vario.data <- vario.calc(db.sub, dirvect=dirvect, opt.code=1, tolcode=0)
# Model using structures that are linear
# at origin (spherical or exponential)
model.vario <- model.auto(vario.data)# , struct=c(1,3,3,2,2))
# Co-Kriging
neigh.kri <- neigh.create(type=2,ndim=3,nmini=nmini,nmaxi=nmaxi,radius=radius,
                          flag.sector=flag.sector,nsect=nsect,nsmax=nsmax,
                          flag.continuous=flag.continuous,flag.aniso=TRUE)

###############
# 4.3 Kriging #
###############
for (survey in unique(db.sub@items$survey_id)){
  message(Sys.time(),": Processing ",survey)
  kri.1 <- kriging(db.sel(db.sub, survey_id==survey), db.grid, model.vario, 
                   neigh.kri, flag.linked=flag.linked)
  beep()
  # Truncating estimations within [0,1]
  ranks = db.ident(kri.1,names="Kriging.Indicator*")
  for(i in ranks){
    kri.1[,i][kri.1[,i] < 0] <- 0
    kri.1[,i][kri.1[,i] > 1] <- 1
  }
  # Building indicators for intervals
  res <- kri.1
  natt.first = res$natt + 1
  res <- db.add(res, zero = (1-res[,ranks[1]]))
  res <- db.add(res, low = (res[,ranks[1]]-res[,ranks[2]]))
  res <- db.add(res, medium = (res[,ranks[2]]-res[,ranks[3]]))
  res <- db.add(res, large = (res[,ranks[3]]-res[,ranks[4]]))
  res <- db.add(res, extrem = (res[,ranks[4]]))
  # Getting rank of most probable interval
  res <- db.compare(res,fun="maxi",names=natt.first:res$natt)
  res <- db.add(res, class=rep(NA,res@nx[1]*res@nx[2]))
  for(i in 1:5){
    res[,"class"][res[,natt.first+i-1]==res[,"maxi"]] <- i
  }
  res <- kres_UTMkm2UTMm(res)
  saveRDS(res, file=paste0(model_out_dir,"3Dresdb_indicator_",survey,".rds"))
  # make dataframe for plotting
  res.df <- as.data.frame(res@items[res@items$Polygon == TRUE,])
  res.df$survey_id=survey
  if (which(unique(db.sub@items$survey_id) == survey) == 1){
    result <- res.df
  } else {
    result <- rbind(result, res.df)
  }
}

################
# 5.1 Plotting #
################
basemap <- load.basemap('montague', zoom=12)
if (make_plots){
  levs <- c("Zero","Low","Medium","High","Extreme")
  result$class <- factor(levs[result$class], levels=levs)
  # save dataframe
  saveRDS(res, file=paste0(model_out_dir,"3Dresdf_indicator_allSurveys.rds"))
  subdat = result[result$x3 %in% c(10,20,30,40),]
  p <- ggplot(subdat, aes(x=x1, y=x2, fill=class)) +  
    geom_raster() +
    scale_fill_manual(values=cmocean(cmap)(5)[1:5], na.translate = F, drop=FALSE) +
    facet_grid(cols=vars(x3), rows=vars(survey_id)) +
    xlab("UTM Eastings (km)")+ylab("UTM Northings (km)")
  ggsave(paste0(plots_out_dir,"3Dindicator_kriging.png"), p, width=10, height=14, dpi=300)
}