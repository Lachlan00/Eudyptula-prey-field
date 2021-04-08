# Result plots draft
source('scripts/eudyptula.R')
source('scripts/visual.R')
source('scripts/HMM/HMM_functions.R')

#####################
# HMM Model results #
#####################
# Load models
mod.hmm <- readRDS('data/simulations/HMM_models/HMM.rds')
hmm.trackSource <- readRDS('data/simulations/HMM_models/HMM_trackSource.rds')

# __Plot model parameters__ #
plot(mod.hmm, breaks=20) # ignore the tracks plots
# __Plot states__ #
plot.HMM.2state(mod.hmm, type='split_points')
# __Plot state desnity __ #
plot.HMM.2state.density.survey(mod.hmm, hmm.trackSource$survey_id, zoom=10)

######################
# Simulation results #
######################
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
tracks.sim <- readRDS('./data/analysis_datasets/cumsum/tracks_sim.rds')
# Filter tracks that leave the survey area like in previous analysis
tracks.real <- inside.survey.zone(tracks.real, threshold=10, 
                                  plot.map=F, plot.title='Real Tracks')
tracks.sim <- inside.survey.zone(tracks.sim, threshold=10,
                                 plot.map=F, plot.title='Simulated Tracks')
ids <- sample(unique(tracks.sim$id), length(unique(tracks.real$id)))
data  <- tracks.sim[tracks.sim$id %in% ids,]
data <- rows.and.levels(data)
# lonlat lims
custom.limits <- list(xmin=min(c(data$lon, tracks.real$lon), na.rm=T),
                      xmax=max(c(data$lon, tracks.real$lon), na.rm=T),
                      ymin=min(c(data$lat, tracks.real$lat), na.rm=T),
                      ymax=max(c(data$lat, tracks.real$lat), na.rm=T))
# Filter tracks to subset
p.real <- tracks.map(tracks.real, legend = F, zoom=10, title='Real Tracks', custom.limits=custom.limits)
p.sim <- tracks.map(data, legend = F, zoom=10, title='Simulated Tracks', custom.limits=custom.limits)
ggarrange(p.real, p.sim)

######################
# Vertcial structure #
######################
# Paramaters
depth = 30 # Depth to extract Sv_mean data
# load databases
surveys <- c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
dbs.Svmn <- load.krig.dbs(var='Sv_mean', surveys=surveys, return_items=T)
# Extract the top water of water column to depth
dbs.Svmn <- lapply(dbs.Svmn, function(db) db[db$x3 <= depth & db$Polygon,])
# Take the mean over x3 (depth)
dbs.Svmn.vert <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=mean))
# Add survey name to each database set
for (i in 1:length(surveys)){
  dbs.Svmn[[i]]$survey <- surveys[i]
  dbs.Svmn.vert[[i]]$survey <- surveys[i]
}
# Plot!!!
ggplot(do.call(rbind, dbs.Svmn.vert), 
                aes(x=x3, y=Kriging.Sv_mean.estim, col=survey)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  coord_flip() +
  labs(x='Depth (m)', y='Sv mean') +
  scale_color_discrete(name='Year', labels=as.character(2015:2019))

################
# Mean density #
################
for (i in 1:length(dbs.Svmn)) {
  message(dbs.Svmn[[i]]$survey[1])
  print(mean(dbs.Svmn[[i]]$Kriging.Sv_mean.estim, na.rm=T))
  print(sd(dbs.Svmn[[i]]$Kriging.Sv_mean.estim, na.rm=T))
}
ggplot(do.call(rbind, dbs.Svmn), aes(x=survey, y=Kriging.Sv_mean.estim))  +
  geom_boxplot(outlier.shape = NA)

################
# Dive density #
################
dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > 5,]

# Mean depth
p.d.mn <- ggplot(dive.data.long, aes(x=depth_mean, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Mean depth (m)') +
  ggtitle('(a)')

# Max depth
p.d.mx <- ggplot(dive.data.long, aes(x=depth_max, col=factor(year(start)))) +
  geom_density() + 
  labs(col='Year') +
  xlab('Max depth (m)') +
  ggtitle('(b)')

# Duration
p.d.dr <- ggplot(dive.data.long, aes(x=dive_duration, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Duration (s)') +
  ggtitle('(c)')

ggarrange(p.d.mn,  p.d.mx, p.d.dr,common.legend=T, legend='right')

############################
# Count raw track fix data #
############################
tracks.raw <- readRDS('data/gps/penguin_tracks_all_CLEAN.rds')
postFilt.ids <- readRDS('data/analysis_datasets/simulations/trackFiltIDs.rds')
tracks.raw <- tracks.raw[tracks.raw$id %in% postFilt.ids,]
tracks.raw <- rows.and.levels(tracks.raw)
tracks.raw <- tracks.time(tracks.raw)
summary(tracks.raw$time)
sd(tracks.raw$time, na.rm=T)
tracks.raw <- rows.and.levels(tracks.raw)

# more stats
tracks.real <- rows.and.levels(tracks.real)
stats.real <- peng.summary.stats(tracks.real)
summary(stats.real)
sd(stats.real$duration)
sd(stats.real$distance_total)
sd(stats.real$displacement)

# report for each year
for (res in split(stats.real, stats.real$survey_id)){
  message(as.character(res$survey_id[1]))
  print(summary(res))
}


