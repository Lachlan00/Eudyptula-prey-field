# Survey area map with penguin tracks
library(rgdal)
library(stringr)
library(ggplot2)
library(ggpubr)
library(maps)
library(sp)
library(maptools)
source('scripts/eudyptula.R')

#####################################################
# Study region, penguin tracks and simulated tracks #
#####################################################
# --- Data Prep --- #
# Load in transect data
transect.lines <- read.csv('data/transects/transect_lines.csv')
# Reshape for ggplot
transect.lines1 <- transect.lines[,1:3]
transect.lines2 <- transect.lines[,c(1,4,5)]
names(transect.lines2)  <- names(transect.lines1)
transect.lines <- rbind(transect.lines1, transect.lines2)
# Load in coastline shapefile
fn = './assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84'
landShp <- readOGR(dsn=path.expand(str_sub(fn,1,-(nchar(basename(fn)))-1)), 
                   layer=basename(fn))
# Load in survey area
fn = "./assets/shp/survey-coast/survey_coast_straight"
surveyShp <- readOGR(dsn=path.expand(str_sub(fn,1,-(nchar(basename(fn)))-1)), 
                layer=basename(fn))
# remove Island
surveyShp@polygons[[1]]@Polygons[[2]] <- NULL
# Load in penguin and simulated tracks
# Load the data
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
tracks.sim <- readRDS('./data/analysis_datasets/cumsum/tracks_sim.rds')
# Filter tracks that leave the survey area like in previous analysis
tracks.real <- inside.survey.zone(tracks.real, threshold=20, use.percentage = T, 
                                  plot.map=F, plot.title='Real Tracks')
tracks.sim <- inside.survey.zone(tracks.sim, threshold=20, use.percentage = T,
                                 plot.map=F, plot.title='Simulated Tracks')
# subset the simulated tracks
# N = length(unique(tracks.real$id))
# ids <- sample(unique(tracks.sim$id), N)
# tracks.sim <- tracks.sim[tracks.sim$id %in% ids,]
# sim.count <- length(unique(tracks.sim$id))

# Make joint dataframe
tracks <- rbind(tracks.real[,c('id','lon','lat','type')], 
                tracks.sim[,c('id','lon','lat','type')])

# Report numbers
message('Real: ',length(unique(tracks.real$id)))
message('Sim: ',sim.count)

# --- Plot Map --- #
xlim <- c(150.04, 150.3)
ylim <- c(-36.47, -36.15)
p.surveyMap <- ggplot() +
  # penguins
  geom_path(data=tracks, mapping=aes(x=lon, y=lat, col=type, group=id), alpha=.15, size=.5,
            show.legend=FALSE) +
  geom_point(data=tracks, mapping=aes(x=lon, y=lat, col=type), alpha=.075, size=.9,
             show.legend=FALSE) +
  scale_color_manual(values=c("#d35400", "#2980b9")) +
  # transects
  geom_line(data=transect.lines, mapping=aes(x=lon1, y=lat1, group=transect), size=.75) +
  # coastline
  geom_polygon(data=landShp, aes(x=long, y=lat, group=group), 
               fill='grey', col='#636363', size=.6) +
  # survey area
  geom_polygon(data=surveyShp, aes(x=long, y=lat, group=group), 
               fill='transparent', col='#c0392b', size=.7) +
  geom_polygon(data=surveyShp, aes(x=long, y=lat, group=group), 
               fill='transparent', col='#c0392b', size=.9, linetype='dotted') +
  # settings
  coord_map(xlim=xlim,  ylim=ylim) +
  theme_bw() + grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL) +
  theme(text = element_text(size=20)) +
  facet_wrap(~type) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# --- Make Map Inset --- #
# map rect extent
zone <- data.frame(x=mean(xlim), y=mean(ylim))
# load world
world <- map("world", fill=TRUE, plot=FALSE)
IDs <- sapply(strsplit(world$names, ":"), function(x) x[1])
world <- map2SpatialPolygons(world, IDs=IDs,
                             proj4string=CRS("+proj=longlat +datum=WGS84"))
world_map <- fortify(world)

p.inset <- ggplot() +
  geom_map(data=world_map, map=world_map,
           mapping = aes(x=long, y=lat, map_id=id), 
           fill="grey", color='#636363') +
  xlim(c(110,160)) + ylim(c(-44,-11)) +
  coord_quickmap() +
  theme_bw() + grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL) +
  theme(text = element_text(size=14)) +
  scale_x_continuous(breaks = pretty(110:160, n = 2), lim=c(110,160)) +
  geom_point(data=zone, mapping=aes(x=x, y=y), col='#c0392b', size=3) +
  theme(panel.grid.minor = element_blank())

# Style
# ggplot transparency
gg_transaprent <- theme(plot.background = element_rect(fill = "transparent", color = NA))
  
# Save files
ggsave(p.surveyMap + gg_transaprent, filename='output/paper_figures/simulation_paper/mapMain_main.png',
       bg="transparent",dpi=300, width=3000/300, height=3000/300)
ggsave(p.inset + gg_transaprent, filename='output/paper_figures/simulation_paper/mapMain_inset.png',
       bg="transparent",dpi=300, width=600/300, height=600/300)

