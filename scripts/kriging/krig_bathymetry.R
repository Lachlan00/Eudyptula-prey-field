# kriging bathymetry
library(RGeostats)
library(cmocean)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(stringi)
library(stringr)
library(ncdf4)
library(ggnewscale)

setwd("/Users/lachlanphillips/Development/PhD/repos/Eudyptula")
source('./scripts/utilities.R')

############
# Settings #
############
input_dir <- './data/bathymetry/'
polygon_fn <- './kriging/assets/kriging_field/survey-coast-straight_poly.rds'
output_dir <- './kriging/output/acoustics/bathymetry/'
xy_res = 3 # rounding factor for lat/lon
convert_UTM = TRUE
gridsize <- c(200, 200) # UTM
# options
save_model <- TRUE # save model as RDS
save_plots <- TRUE # save plots?
animate <- TRUE # make animation?
colour_grid <- TRUE # color grid?

###############################
# Add Hydrodynamic Bathy data #
###############################
# leave this FALSE as it doesn't work
add_hydro_bathy <- FALSE
ncfname <- '/Volumes/LP_MstrData/master-data/ocean/ROMS/highres/sample/out_his_01551.nc'

#############
# Load data #
#############
# read in 2016 bathymetry data
seabed.file.ls <- list.files(path=input_dir, pattern=".csv$", recursive=FALSE)
df.ls <- list()
i <- 1
for (file in paste0(input_dir, seabed.file.ls)){
  df.ls[[i]] <- read.csv(file)
  i <- i + 1
}
seabed.df <- do.call('rbind', df.ls)
# keep only data we want
seabed.df <- seabed.df[,c('Longitude','Latitude','Depth')]
colnames(seabed.df) <- c('lon', 'lat', 'depth')
seabed.df$source <- 'acoustics'

# also read in the hydrodynamic model data
if (add_hydro_bathy){
  ncin <- nc_open(ncfname)
  nc_bathy <- as.vector(ncvar_get(ncin, 'h'))
  nc_lons <- as.vector(ncvar_get(ncin, 'lon_rho'))
  nc_lats <- as.vector(ncvar_get(ncin, 'lat_rho'))
  # make data frame
  nc_df <- data.frame(lon=nc_lons, lat=nc_lats, depth=nc_bathy)
  nc_df$source <- 'ROMS'
  # add to seabed dataframe
  seabed.df <- rbind(seabed.df, nc_df)
}

# add custom points
custom_depths <- read.csv('./data/bathymetry/custom_depths/baranguba_depths.csv')
custom_depths$source <- 'custom'
seabed.df <- rbind(seabed.df, custom_depths)

# take the mean by location to reduce processing time
seabed.df <- aggregate(seabed.df, by=list(round(seabed.df$lon, xy_res), 
                                          round(seabed.df$lat, xy_res)),
                       FUN="mean", na.action=na.pass)
# organise
seabed.df <- seabed.df[,c(1:2,5)]
colnames(seabed.df) <- c('lon', 'lat', 'depth')

###########
# Kriging #
###########
# get lon and lat in UTM
if (convert_UTM){
  seabed.df <- lonlat2UTM.df(seabed.df)
}
# make database
db.floor <- db.create(seabed.df)
db.floor <- db.locate(db.floor, c("lon","lat"),"x")
db.floor <- db.locate(db.floor, 'depth',"z")

# read polygon
poly <- readRDS(polygon_fn)
# Correct polygon Lon and Lat to be in UTM
if (convert_UTM){
  poly <- lonlat2UTM.poly(poly)
}
#initialise a grid
db.grid <- db.grid.init(poly, dcell=gridsize) # nodes=c(100,50,20))
# make grid from polygon
# db.grid <- db.polygon(db.grid, poly)
#Calculate experimental variogram
vari <- vario.calc(db.floor)
# Fit variogram model
vg.mod <- model.auto(vari, draw=TRUE)
# Define neighbourhood
neimov <- neigh.create(type=0,ndim=2) # nmini=7,nmaxi=15, radius=10000
# Run the kriging
kres <- kriging(dbin=db.floor, dbout=db.grid, model=vg.mod, neigh=neimov)
plot(kres)

if (save_model){
  saveRDS(kres, file=paste0(output_dir,"models/krig_bathymetry-db.rds"))
}
# Convert kriging results into dataframe
t.df <- as.data.frame(kres@items)
names(t.df)<- c('rank', 'lon','lat','polygon','depth','stdev')
# Convert Lon lat into degrees
if (convert_UTM){
  t.df <- UTM2lonlat.df(t.df)
}
#Set all data outside of Polygon to be NA
t.df[,'depth'][t.df$polygon==FALSE] <- NA
if (save_model){
  saveRDS(t.df, file=paste0(output_dir,"models/krig_bathymetry-df.rds"))
}

############
# Plotting #
############
# geom_raster tries to create a giant plot because the resolution of the plot is based on the 
# minimal distance between points. Consequently it gets a vector memory error. To avoid this
# we plot kres@items with UTM coordinates and then reporject lat and lon on after.
# NOTE: For now make plot in UTM coordinates
# make a UTM dataframe
UTM_df <- kres@items
UTM_df[,'Kriging.depth.estim'][UTM_df$Polygon==FALSE] <- NA
names(UTM_df) <- c('rank','lon','lat','poly','depth','stdev')

depth_plot <- ggplot(data=UTM_df, aes(x=lon, y=lat, fill=depth)) +
  geom_raster() +
  scale_fill_gradientn(colours=cmocean('deep')(100)) +
  labs(x='Easting', y='Northing',fill='Depth') +
  ggtitle('Kriged Bathymetry - UTM Projection')

if (save_plots){
  ggsave(paste0(output_dir,"plots/krig_bathymetry.png"),
         depth_plot, width=7, height=7, dpi=200)
}

depth_plot <- depth_plot +
  geom_point(data=seabed.df, size=2, pch=21)

if (save_plots){
  ggsave(paste0(output_dir,"plots/krig_bathymetry_datasource.png"),
         depth_plot, width=7, height=7, dpi=200)
}

# 3d look
# make dpeth into matrix
depth_matrix <- matrix(t.df$depth, nrow=kres@nx[1], ncol=kres@nx[2])*-1
inv_depth_matrix <- depth_matrix[,c(kres@nx[2]:1), drop=FALSE]
Depth3D <- plot_ly(showscale=FALSE) %>%
  layout(
    autorange = F, 
    aspectmode = 'manual', 
    scene = list(
      aspectratio = list(x = 1, y = 1, z = 0.35)
    )
  ) %>%
  add_surface(z = inv_depth_matrix, colors=cmocean('deep')(100)[100:1])
# SAVE
if (save_plots){
  saveWidgetFix(Depth3D, paste0(output_dir,"plots/krig_bathymetry_3D.html"), selfcontained=TRUE)
}

# perspective plot
# make colours
depth_round <- round(depth_matrix,0)*-1
cols_set <- rev(cmocean('deep')((max(depth_round,na.rm=T)+1) - min(depth_round,na.rm=T)))
zfacet <- depth_matrix[-1, -1] + depth_matrix[-1, -kres@nx[2]] + 
  depth_matrix[-kres@nx[1], -1] + depth_matrix[-kres@nx[1], -kres@nx[2]]
facetcol <- cut(zfacet, (max(depth_round,na.rm=T)+1) - min(depth_round,na.rm=T))
depth_cols <- cols_set[facetcol]

if (!colour_grid){
  depth_cols <- 'white'
}

# grid plot
if (save_plots){
  png(paste0(output_dir,"plots/krig_bathymetry_persp.png"),
      width=1000, height=850)
}
persp(depth_matrix, theta=55, phi=25, expand=0.3, xlab='lon', ylab='lat', 
        zlab='Depth', col=depth_cols)
if (save_plots){
  dev.off()
}

# animation
if (animate){
  # make rotating animation
  rand_id <- stri_rand_strings(1, 8)
  dir.create(paste0(output_dir,'/plots/anim_',rand_id))
  # make pngs
  for (i in  1:359){
    png(paste0(output_dir,"plots/anim_",rand_id,"/krig_bathymetry_persp_",str_pad(i, 3, pad = "0"),".png"),
        width=750, height=500)
    persp(depth_matrix, theta = i-1, phi = 23, expand = 0.35, xlab='lon', ylab='lat', 
            zlab='depth', col=depth_cols)
    dev.off()
    message(paste('Processing frame',i,'of 359..'), "\r", appendLF=FALSE)
    flush.console()
  }
  message('')
  message('Joining frames..')
  # merge
  abs_path <- normalizePath(paste0(output_dir,'/plots'))
  system(paste0('convert ',abs_path,'/anim_',rand_id,'/*.png ',abs_path,'/bathy.gif'))
  message('Converting frame delay..')
  system(paste0('convert -delay 5x100 ',abs_path,'/bathy.gif ',abs_path,'/bathy.gif'))
  # delete components
  message('Cleaning up..')
  unlink(paste0(output_dir,'plots/anim_',rand_id,'/'), recursive=TRUE)
}

