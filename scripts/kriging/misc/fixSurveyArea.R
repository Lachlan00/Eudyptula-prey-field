# Fix up kriging field
# This script will produce a shapefile with variable dimensions based on the transect area 
library(sp)

# Load in transect lines
transect <- read.csv('data/transects/transect_lines.csv')
# Reshape for ggplot
transect.lines1 <- transect[,1:3]
transect.lines2 <- transect[,c(1,4,5)]
names(transect.lines2)  <- names(transect.lines1)
transect.lines <- rbind(transect.lines1, transect.lines2)

# Load in coastline shapefile
fn = './assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84'
landShp <- readOGR(dsn=path.expand(str_sub(fn,1,-(nchar(basename(fn)))-1)), 
                   layer=basename(fn))

# Calulate difference between latitudes (it's consistent)
# We will add half the value to the top and bottom and east of the range
lat.diff <-  min(diff(transect$lat1))
val <- lat.diff/2

# Setup north points  (west to east)
N.lon <- c(150.05, transect$lon2[nrow(transect)])
N.lat <- rep(max(transect$lat1) + val, 2)
# Setup east points (north to south)
E.lon <- rev(transect$lon2) + add
E.lat <- rev(transect$lat2)
# Setup South points (east to west)
S.lon <- c(transect$lon2[1], 150.05)
S.lat <- rep(min(transect$lat1) - val, 2)
# Make dataframe
shpPoints <- data.frame(lon=c(N.lon, E.lon, S.lon),
                        lat=c(N.lat, E.lat, S.lat))

# Check all is A-okay
ggplot() + 
  geom_line(aes(x=lon1, y=lat1, group=transect), transect.lines, size=.75) +
  geom_path(aes(x=lon, y=lat), shpPoints, color='red') +
  geom_point(aes(x=lon, y=lat), shpPoints, color='red') +
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  coord_map(xlim=c(150.04, 150.3),  ylim=c(-36.47, -36.15)) +
  theme_bw() + grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL) +
  theme(text = element_text(size=14))

# Make a shape file from the points
names(shpPoints) <- c('long','lat')
shapeOut <- Polygon(shpPoints)
shapeOut <- Polygons(list(shapeOut),1)
shapeOut <- SpatialPolygons(list(shapeOut))
proj4string(shapeOut)<- CRS("+proj=longlat +datum=WGS84")
raster::shapefile(shapeOut, "assets/shp/kriging_surveyArea/kriging_surveyArea.shp")

# Difference with coast calculated in QGIS 3

