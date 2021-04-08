# Kriging area and HMM overlay
# libraries
library(pbapply)
library(plot3D)
library(plot3Drgl)
library(ggplot2)
library(cmocean)
library(beepr)
library(tools)
library(raster)
library(momentuHMM)

source('scripts/Eudyptula.R')

weight.depth = TRUE

# If using a weighted mean load dive data
if (weight.depth){
  # load penguin dive data
  dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
  # Load real tracks to filter with
  tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
  # dives longer than 5 seconds
  longCutoff <- 5
  dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > longCutoff,]
  dive.data.long$year <- year(dive.data.long$start)
  dive.data.long <- dive.data.long[as.character(dive.data.long$id) %in% 
                                     as.character(tracks.real$id),]
  dive.data.long <- tracks.nearest.survey(dive.data.long, dtUTCcol='start')
  dive.data.long <- rows.and.levels(dive.data.long, idcol='survey_id')
  # Split and calculate the distribution
  dive.dist <- data.frame(survey_id=dive.data.long$survey_id,
                          depth=as.integer(round(dive.data.long$depth_max)),
                          duration=dive.data.long$dive_duration)
  dive.dist.ls <- split(dive.dist, dive.dist$survey_id)
  for (i in 1:length(dive.dist.ls)){
    res <- data.frame(depth=0:max(dive.dist$depth),
                      freq=0)
    tab <- table(dive.dist.ls[[i]]$depth)
    res[as.integer(names(tab))+1,'freq'] <- tab
    res$survey_id <- dive.dist.ls[[i]]$survey_id[1]
    dive.dist.ls[[i]] <- res
  }
  dive.weights <- do.call(rbind, dive.dist.ls)
  weight.df <- rows.and.levels(dive.weights, idcol='survey_id')
}

###########
# 3D plot #
###########
# cmap
cmap <- rev(cmocean('deep')(100))[1:100]
# Load all databases
surveys = c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
krig.dfs <- load.krig.dbs('Sv_mean', surveys, return_items = T)
# Reduce all to the top 30m or top weight level
if (weight.depth){
  krig.dfs <- lapply(krig.dfs, function(db) db[db$x3 <= max(weight.df$depth),])
} else {
  krig.dfs <- lapply(krig.dfs, function(df) df[df$x3 <= 30,])
}

# outside Polygon gets NA
for (i in 1:length(krig.dfs)){
  krig.dfs[[i]][!krig.dfs[[i]]$Polygon,'Kriging.Sv_mean.estim'] <- NA
}

# load a single database for 3d plot
krig.df <- krig.dfs[[2]]

# Make 2d
if (weight.depth){
  # Weight by dive count
  for (i in 1:length(krig.dfs)){
    krig.dfs[[i]] <- aggregate(krig.dfs[[i]], by=list(krig.dfs[[i]]$x1, krig.dfs[[i]]$x2),
                          FUN = function(x) weighted.mean(x, weight.df$freq[weight.df$survey_id == names(krig.dfs)[i]]))
    krig.dfs[[i]]$survey <- surveys[i]
  }
} else {
  # Or make all dfs a the mean to 30m
  for (i in 1:length(krig.dfs)){
    krig.dfs[[i]] <- aggregate(krig.dfs[[i]], 
                               by=list(krig.dfs[[i]]$x1, krig.dfs[[i]]$x2), 
                               FUN=mean)
    krig.dfs[[i]]$survey <- surveys[i]
  }
}


krig.dfs <- do.call(rbind, krig.dfs)
names(krig.dfs)[1:2] <- c('x','y')
# convert x and y to lon/lat (approximations)
krig.dfs$lon <- UTM2lonlat(krig.dfs$x, rep(krig.dfs$y[1], length(krig.dfs$x)))$x
krig.dfs$lat <- UTM2lonlat(rep(krig.dfs$x[1], length(krig.dfs$y)), krig.dfs$y)$y

# Find min and max Sv
vmin <- ceiling(min(c(min(krig.dfs$Kriging.Sv_mean.estim, na.rm=T), 
                      min(krig.df$Kriging.Sv_mean.estim, na.rm=T))))
vmax <- floor(max(c(max(krig.dfs$Kriging.Sv_mean.estim, na.rm=T),
                    max(krig.df$Kriging.Sv_mean.estim, na.rm=T))))

# manual vmin and vmax
vmin <- -67
vmax <- -38

krig.dfs$survey <- as.factor(substr(krig.dfs$survey,1,4))

#Set up dimension
x <- unique(krig.df$x1)
y <- unique(krig.df$x2)
# convert x and y to lon/lat (approximations)
lon <- UTM2lonlat(x, rep(y[1], length(x)))$x
lat <- UTM2lonlat(rep(x[1], length(y)), y)$y
z <- -unique(krig.df$x3)
#Create 3D data array
p <- array(krig.df$Kriging.Sv_mean.estim, dim=c(length(x),length(y),length(z)))

# Plot
png('output/paper_figures/simulation_paper/krigMap_3d.png',
    width=800, height=800)
p.3d <- slice3D(lon, lat, z, colvar = p,
        col = cmap,
        colkey=F,
        NAcol = "transparent",
        clim = c(vmin, vmax),
        xs = NULL,
        ys = lat[c(1)],
        zs = c(z[c(1,8,18,30)]),
        d = 1,
        xlab="Longitude",
        ylab="Latitude",
        zlab="Depth (m)",
        clab='Sv mean',
        ticktype="detailed",
        cex.axis=1.3,
        cex.lab=1.3)
dev.off()

#############################
# Plot the maps of mean 30m #
#############################
# load HMM models
mod.HMM <- readRDS('./data/simulations/HMM_models/HMM.rds')
tracks <- mod.HMM$data
tracks$state <- viterbi(mod.HMM)
names(tracks)[c(1,4,5,6)] <- c('id','lon','lat','survey')
# Filter tracks that leave the survey area like in previous analysis
tracks <- inside.survey.zone(tracks, threshold=10, 
                                  plot.map=F, plot.title='Real Tracks')
# get only foraging state
tracks <- tracks[tracks$state == 1,]
tracks$survey <- as.factor(substr(tracks$survey,1,4))

# Load in coastline shapefile
fn = './assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84'
landShp <- readOGR(dsn=path.expand(str_sub(fn,1,-(nchar(basename(fn)))-1)), 
                   layer=basename(fn))
# Map limits
xlim <- c(min(krig.dfs$lon), max(krig.dfs$lon))
ylim <- c(min(krig.dfs$lat), max(krig.dfs$lat))

# Transform polygon
# Hacky method
# convert to UTM and then back to lon lat using the same approximation as before
landShpT <- spTransform(landShp, CRS("+proj=utm +zone=56H +ellps=WGS84"))
landShp.Store <- landShp
for (i in 1:length(landShpT@polygons)){
  for (j in length(landShpT@polygons[[i]]@Polygons)){
    # Convert lon
    landShp@polygons[[i]]@Polygons[[j]]@coords[,1] <-  
      UTM2lonlat(landShpT@polygons[[i]]@Polygons[[j]]@coords[,1], 
                 rep(krig.dfs$y[1], 
                     length(landShpT@polygons[[i]]@Polygons[[j]]@coords[,1])))$x
    # Convert lat
    landShp@polygons[[i]]@Polygons[[j]]@coords[,2] <-  
      UTM2lonlat(rep(krig.dfs$x[1], 
                     length(landShpT@polygons[[i]]@Polygons[[j]]@coords[,2])),
                 landShpT@polygons[[i]]@Polygons[[j]]@coords[,2],)$y
  }
}

# # Plot
# p.maps <- ggplot(krig.dfs, aes(x=lon, y=lat, fill=Kriging.Sv_mean.estim)) + 
#   geom_raster() +
#   scale_fill_gradientn(colours=cmap, limits=c(vmin, vmax), na.value="grey66",
#                        guide=guide_colorbar(barwidth=25, frame.colour=c("black"))) +
#   facet_wrap(~survey) +
#   theme_bw() + grids(linetype = "dashed") +
#   labs(x=NULL, y=NULL, alpha=NULL, size=NULL, fill=NULL) +
#   theme(text = element_text(size=14)) +
#   theme(legend.position = "bottom") +
#   coord_cartesian(xlim=xlim,  ylim=ylim, expand=c(0,0)) +
#   stat_density_2d(data=tracks, mapping=aes(x=lon, y=lat), col='white', alpha=.5,
#                   inherit.aes = F, geom = "density2d") +
#   # coastline
#   geom_polygon(data=landShp, mapping=aes(x=long, y=lat, group=group),
#                fill='grey', col='#636363', size=.6, alpha=1, inherit.aes = F)

coast <- raster::crop(landShp, extent(xlim, ylim))

p.maps <- ggplot()  +  
  geom_raster(data = krig.dfs, aes(x=lon, y=lat, fill=Kriging.Sv_mean.estim)) +
  scale_fill_gradientn(colours=cmap, limits=c(vmin, vmax), na.value="grey66",
                       guide=guide_colorbar(barwidth=25, frame.colour=c("black"))) +
  theme_bw() + #grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL, fill=NULL) +
  theme(text = element_text(size=17)) +
  theme(legend.position = "bottom") +
  geom_density_2d(data=tracks, mapping=aes(x=lon, y=lat), col='white', alpha=.5,
                  inherit.aes = F) +
  geom_polygon(data=coast, mapping=aes(x=long, y=lat, group=group),
               fill='grey', col='#636363', size=.6, alpha=1, inherit.aes = F)+
  facet_wrap(~survey) +
  coord_cartesian(xlim=xlim,  ylim=ylim, expand=c(0,0))

# Save p.maps
ggsave(p.maps, filename='output/paper_figures/simulation_paper/krigMap_pengs.png',
       bg="transparent",dpi=300, width=3500/300, height=3000/300)


