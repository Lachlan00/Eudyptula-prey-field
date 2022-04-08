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
library(directlabels)

source('scripts/Eudyptula.R')

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
xlim <- c(min(tracks$lon, na.rm=T)-.2, max(tracks$lon, na.rm=T)+.2)
ylim <- c(min(tracks$lat, na.rm=T)-.2, max(tracks$lat, na.rm=T)+.2)

coast <- raster::crop(landShp, extent(xlim, ylim))

# remove island
coast <- subset(coast, area > 10)

# make spatial points
points <- as.data.frame(tracks[,c('lon','lat')])
points <- points[complete.cases(points),]
coordinates(points) <- ~ lon + lat

# Fix coords
proj4string(points) <- CRS("+proj=longlat +datum=WGS84")
proj4string(coast) <- CRS("+proj=longlat +datum=WGS84")

plot(coast)
plot(points, add=T)

# calc distance
library(rgeos)
dists <- apply(gDistance(points, coast, byid=T), 2, min)
# approximate 1 degree longitude
#mean(tracks$lat, na.rm=T) # = -36.34521 therefore lon = ~ 89 kilometers 
dists <- as.vector(dists*89)

# attach to tracks
tracks.dists <- tracks[!is.na(tracks$lon),]
tracks.dists$dists <- dists

# calc mean and sd for each year
for (df in split(tracks.dists, tracks.dists$survey)){
  message('\n',as.character(df$survey[1]),' | mean: ', round(mean(df$dists),2),
          ' ± ',round(sd(df$dists),2),' km')
  message(paste('quantile:',round(quantile(df$dists,probs=c(0.025,0.975)),2)))
}

# plot 
p.dense <- ggplot(tracks.dists, aes(x=dists, fill=survey, color=survey)) +
  geom_density(alpha=0.3) + 
  labs(x='Distance from coast (km)', y='Density', fill='Year', color='Year')
  
