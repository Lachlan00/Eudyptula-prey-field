# Post hoc analysis to see why 2017 is different
library(ggplot2)
library(cmocean)
library(scales)
library(sp)
library(momentuHMM)
library(beepr)
source('scripts/eudyptula.R')
source('scripts/kde_percent.R')

# Paramaters
depth = 30 # Depth to extract Sv_mean data
kde.per = 75

# load databases
surveys <- c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
dbs.Svmn <- load.krig.dbs(var='Sv_mean', surveys=surveys, return_items=T)
# Extract the top water of water column to depth
dbs.Svmn <- lapply(dbs.Svmn, function(db) db[db$x3 <= depth & db$Polygon,])

# Side caclulation (report total Sv mean)
for (i in 1:length(dbs.Svmn)){
  message(names(dbs.Svmn)[i])
  print(mean(dbs.Svmn[[i]]$Kriging.Sv_mean.estim))
  print(sd(dbs.Svmn[[i]]$Kriging.Sv_mean.estim))
}

# Aggregate in 2 dimesnions
dbs.Svmn.vert <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=mean))
dbs.Svmn.vert.sd <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=sd))

# Side calcs 
# Peak mean and depth for 2015
dbs.Svmn.vert$`2015_S1`[which.max(dbs.Svmn.vert$`2015_S1`$Kriging.Sv_mean.estim),]
# 2016 8-24 m depth range
mean(dbs.Svmn.vert$`2016_S2`[dbs.Svmn.vert$`2016_S2`$x3 > 7 & 
                             dbs.Svmn.vert$`2016_S2`$x3 <= 23,
                             'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2016_S2`[dbs.Svmn.vert$`2016_S2`$x3 > 7 & 
                               dbs.Svmn.vert$`2016_S2`$x3 <= 23,
                             'Kriging.Sv_mean.estim'])
# 2017 summaries
mean(dbs.Svmn.vert$`2017_S2`[dbs.Svmn.vert$`2017_S2`$x3 > 20, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2017_S2`[dbs.Svmn.vert$`2017_S2`$x3 > 20, 'Kriging.Sv_mean.estim'])
# 2018 summaries
mean(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 <= 5, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 <= 5, 'Kriging.Sv_mean.estim'])
mean(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 >= 20, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 >= 20, 'Kriging.Sv_mean.estim'])
dbs.Svmn.vert$`2018_S2`[which.min(dbs.Svmn.vert$`2018_S2`$Kriging.Sv_mean.estim),]
dbs.Svmn.vert$`2018_S2`[which.max(dbs.Svmn.vert$`2018_S2`$Kriging.Sv_mean.estim),]
# 2019
mean(dbs.Svmn.vert$`2019_S2`[dbs.Svmn.vert$`2019_S2`$x3 <= 4, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2019_S2`[dbs.Svmn.vert$`2019_S2`$x3 <= 4, 'Kriging.Sv_mean.estim'])

# Add survey name to each database set
for (i in 1:length(surveys)){
  dbs.Svmn.vert[[i]]$survey <- surveys[i]
  dbs.Svmn.vert[[i]]$sd <- dbs.Svmn.vert.sd[[i]]$Kriging.Sv_mean.estim
}

# load penguin dive data
dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
# dives longer than 5 seconds
longCutoff <- 5
dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > longCutoff,]
dive.data.long$year <- year(dive.data.long$start)

# Load real tracks to filter with
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
dive.data.long <- dive.data.long[as.character(dive.data.long$id) %in% 
                                 as.character(tracks.real$id),]

# Percentage dives < 30 m adn < 10 m
nrow(dive.data.long[dive.data.long$depth_max < 30,])/ nrow(dive.data.long)
nrow(dive.data.long[dive.data.long$depth_max < 10,])/ nrow(dive.data.long)

#####################################
# Plot dive density against Sv mean #
#####################################
# grey col
colGrey <- '#757575'
# Plot side by side
Sv_vert <- do.call(rbind, dbs.Svmn.vert)
Sv_vert$year <- as.numeric(substr(Sv_vert$survey,1,4))

#--- SAME AXIS ---#
# try combine
p.vert <- ggplot(data=NULL) +
  geom_density(data=dive.data.long, mapping=aes(x=depth_max), 
               size=1, col=colGrey, alpha=.34, fill='grey') +
  geom_line(data=Sv_vert, mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17)))), 
            size=1, col='#c0392b', alpha=.8) +
  geom_point(data=Sv_vert, mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17)))), 
            size=2, col='#c0392b', alpha=.8) +
  ylab('Density of dive depths') + xlab('Depth (m)') +
  scale_y_continuous(sec.axis = sec_axis(~rescale(.,
      to=c(min(Sv_vert$Kriging.Sv_mean.estim), max(Sv_vert$Kriging.Sv_mean.estim))), 
      name="Acoustic density (Sv mean)")) +
  coord_flip() +
  scale_x_reverse(limits=c(30,0)) +
  facet_wrap(~year, strip.position='right') +
  theme_bw() + grids(linetype = "dashed") +
  theme(text = element_text(size=16)) +
  theme(strip.background = element_rect(fill="#f2f2f2", colour="black"),
        strip.placement = "outside", 
        strip.text = element_text(size=16))

ggsave(p.vert, filename='output/paper_figures/simulation_paper/diveVert_sameAxis.png',
       dpi=300, width=3500/300, height=3000/300)


##### ---- Revision Update ---- #####
# Add Diving and Sv mean in kernel density areas
# 1 calculate kernel density region
# load HMM models
mod.HMM <- readRDS('./data/simulations/HMM_models/HMM.rds')
m.tracks <- mod.HMM$data
m.tracks$state <- viterbi(mod.HMM)
names(m.tracks)[c(1,4,5,6,7)] <- c('id','lon','lat','survey','state')
# Filter tracks that leave the survey area like in previous analysis
m.tracks <- inside.survey.zone(m.tracks, threshold=20, use.percentage = T, 
                               plot.map=F, plot.title='Real Tracks')
# get only foraging state
m.tracks <- m.tracks[m.tracks$state == 1,]
m.tracks$survey <- as.factor(substr(m.tracks$survey,1,4))
# Calculate kernel density #####################--------------------------------
# p.kernel <- ggplot()  +  
#   facet_wrap(~survey) +
#   stat_density_2d(data=m.tracks, mapping=aes(x=lon, y=lat),
#                   inherit.aes = F, size=.4, bins=4, contour_var = 'ndensity')  +
#   xlim(149.5, 150.5) + ylim(-36.7, -35.7)
# kernel_data <- ggplot_build(p.kernel)$data[[1]]
# kernel_data$survey <- c(2015:2019)[kernel_data$PANEL]

# NEW
kde.df <- KDE_percentage(m.tracks[complete.cases(m.tracks),c('lon','lat')], 
                         levels=c(kde.per), facet=tracks[complete.cases(m.tracks),c('survey')])
kde.df$survey <- kde.df$facet



# find what levels are for each survey
kernel_data.ls <- split(kde.df, kde.df$survey)
# For each kernel data set extract the outside border and create 
# shape files from them
kernel_poly <- list()
for (i in 1:length(kernel_data.ls)){
  kdata <- kernel_data.ls[[i]]
  # Split and make into shape files
  kdata.ls <- split(kdata, kdata$group)
  # make polygons
  for (j in 1:length(kdata.ls)){
    kdata.ls[[j]] <- Polygon(cbind(kdata.ls[[j]]$x, kdata.ls[[j]]$y))
  }
  kernel_poly[[i]] <- SpatialPolygons(list(Polygons(kdata.ls, 1)))
}
names(kernel_poly) <- 2015:2019

# 2 Extract diving behaviour in area
# get times of points inside zone
tracks.real <- tracks.real[complete.cases(tracks.real),]
tracks.real$survey <- substr(tracks.real$survey_id,1,4)
tracks.real.ls <- split(tracks.real, tracks.real$survey)
for (i in 1:length(tracks.real.ls)){
  poly <- kernel_poly[[i]]
  # check if over
  df <- tracks.real.ls[[i]]
  dat <- data.frame(lon = df$lon,
                    lat = df$lat)
  coordinates(dat) <- ~ lon + lat
  # check if points on ocean
  ocean <- over(dat, poly)
  
  print(ggplot() + 
    geom_polygon(data=poly, aes(x=long, y=lat, group=group)) + 
    coord_fixed(1) +
    geom_point(data=df, aes(x=lon, y=lat, group=ocean, color=ocean), 
               size=2, shape=18, na.rm=TRUE) +
    ggtitle(c(2015:2019)[i]))
  
  # report percentage in kernel region
  message('\n',surveys[i])
  print(round(sum(ocean, na.rm = T)/length(ocean)*100,2))
  
  tracks.real.ls[[i]]$over <- !is.na(ocean)
}
tracks.real <- rows.and.levels(do.call(rbind, tracks.real.ls))
# attach dives to include by finding track and nearest timestamp
dive.data.long.ls <- split(dive.data.long, as.character(dive.data.long$id))
for (i in 1:length(dive.data.long.ls)){
  id <- names(dive.data.long.ls)[i]
  dives <- dive.data.long.ls[[i]]
  track <- tracks.real[tracks.real$id == id,]
  if (nrow(track) == 0){
    message("No GPS for ",id)
    next
  }
  # crude method, kd trees would be way faster
  dives$over <- FALSE
  for (j in 1:nrow(dives)){
    mean_dive_t <- mean(c(dives$start[j], dives$end[j]))
    # nearest track timestamp
    dives$over[j] <- track$over[which.min(abs(track$dtUTC - mean_dive_t))]
  }
  dive.data.long.ls[[i]] <- dives
}
kernel_dives <- do.call(rbind, dive.data.long.ls)
# drop non included dives
kernel_dives <- rows.and.levels(kernel_dives[kernel_dives$over,])

# 3 Extract Sv mean in area
# load databases
surveys <- c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
dbs.Svmn <- load.krig.dbs(var='Sv_mean', surveys=surveys, return_items=T)
# Extract the top water of water column to depth
dbs.Svmn <- lapply(dbs.Svmn, function(db) db[db$x3 <= depth & db$Polygon,])
# For each db filter with kernel polygons
for (i in 1:length(dbs.Svmn)){
  db <- dbs.Svmn[[i]]
  names(db)[2:3] <- c('lon','lat')
  # convert db to latlon
  db <- UTM2lonlat.df(db)
  # Check if over polygons
  poly <- kernel_poly[[i]]
  # check if over
  dat <- data.frame(lon = db$lon,
                    lat = db$lat)
  coordinates(dat) <- ~ lon + lat
  # check if points on kernel
  ocean <- over(dat, poly)
  dbs.Svmn[[i]] <- db[!is.na(ocean),]
}
# now do normal aggregations
dbs.Svmn.vert.kernel <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=mean))
dbs.Svmn.vert.kernel.sd <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=sd))
# Add survey name to each database set
for (i in 1:length(surveys)){
  dbs.Svmn.vert.kernel[[i]]$survey <- surveys[i]
  dbs.Svmn.vert.kernel[[i]]$sd <- dbs.Svmn.vert.kernel.sd[[i]]$Kriging.Sv_mean.estim
}
Sv_vert_kernel <- do.call(rbind, dbs.Svmn.vert.kernel)
Sv_vert_kernel$year <- as.numeric(substr(Sv_vert_kernel$survey,1,4))




# Precalulate min/max for kernel and all Sv
Svmin <- min(c(Sv_vert$Kriging.Sv_mean.estim - Sv_vert$sd,
               Sv_vert_kernel$Kriging.Sv_mean.estim - Sv_vert_kernel$sd))
Svmax <- max(c(Sv_vert$Kriging.Sv_mean.estim + Sv_vert$sd,
               Sv_vert_kernel$Kriging.Sv_mean.estim + Sv_vert_kernel$sd))
Svlim <- c(Svmin, Svmax)           

# PLOT
p.vert.mod <- ggplot(data=NULL) +
  # Diving dist
  geom_density(data=dive.data.long, mapping=aes(x=depth_max), 
               size=1, col=colGrey, alpha=.54, fill='grey') +
  geom_density(data=kernel_dives, mapping=aes(x=depth_max), 
               size=1, col='darkred', alpha=.14, fill='red') +
  # Kernel Sv
  geom_errorbar(data=Sv_vert_kernel, 
                mapping=aes(x=x3, 
                            ymin=rescale(Kriging.Sv_mean.estim-sd, to=c(0, 0.17), from=Svlim),
                            ymax=rescale(Kriging.Sv_mean.estim+sd, to=c(0, 0.17), from=Svlim)),
                width=.5, position=position_dodge(.9), col='darkred', alpha=.3, size=1) +
  geom_line(data=Sv_vert_kernel, 
            mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17), from=Svlim))), 
            size=1, col='darkred', alpha=.8) +
  geom_point(data=Sv_vert_kernel, 
             mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17), from=Svlim))), 
             size=2, col='darkred', alpha=.8) +
  # All Sv
  geom_errorbar(data=Sv_vert, 
                mapping=aes(x=x3, 
                            ymin=rescale(Kriging.Sv_mean.estim-sd, to=c(0, 0.17), from=Svlim),
                            ymax=rescale(Kriging.Sv_mean.estim+sd, to=c(0, 0.17), from=Svlim)),
                width=.5, position=position_dodge(.9), col=colGrey, alpha=.3, size=1) +
  geom_line(data=Sv_vert, 
            mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17), from=Svlim))), 
            size=1, col=colGrey, alpha=.8) +
  geom_point(data=Sv_vert, 
             mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17), from=Svlim))), 
             size=2, col=colGrey, alpha=.8) +
  # Plot settings
  ylab('Density of dive depths') + xlab('Depth (m)') +
  scale_y_continuous(sec.axis = sec_axis(~rescale(.,
                                                  to=Svlim), 
                                         name="Acoustic density (Sv mean)")) +
  coord_flip() +
  scale_x_reverse(limits=c(30,0)) +
  facet_wrap(~year, strip.position='right') +
  theme_bw() + grids(linetype = "dashed") +
  theme(text = element_text(size=16)) +
  theme(strip.background = element_rect(fill="#f2f2f2", colour="black"),
        strip.placement = "outside",
        strip.text = element_text(size=16))






ggsave(p.vert.mod, filename='output/paper_figures/simulation_paper/diveVert_sameAxis_kernel.png',
      dpi=300, width=3500/300, height=3000/300)
beep()

