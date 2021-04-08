source("./scripts/HMM/HMM_functions.R")
source("./scripts/visual.R")
source("./scripts/quality_control.R")
library(beepr)
library(ggcorrplot)
library(cmocean)

# SET
survey_year=2018

#################
# Preprocessing #
#################
# model_output
model_out <- paste0('./data/analysis_datasets/HMM/model_HMM2',survey_year,'.rds')
# Data source
# load processed track data
tracks <- readRDS('./data/analysis_datasets/HMM/tracks_reg5m_processed.rds')
# clean gaps
tracks <- tracks.reg.filter.gaps(tracks)
# Final quality control checks
tracks <- quality.GPS.fix.count(tracks, min_threshold=50, filter=T)

# filter by years (optional)
tracks <- track.date.filter(tracks, paste0(survey_year,'-01-01'), paste0(survey_year+1,'-01-01'))
tracks <- tracks.survey.time.filter(tracks, window=5)

# replace dive statistics NAs with zeros
dive_cols <- c('diveDepth_mean','diveDepth_max','diveDuration_mean',
               'diveDuration_max','diveDuration_sum')
tracks$diveCount[is.na(tracks$diveCount)] <- 0
tracks[is.na(tracks$diveDuration_sum),dive_cols] <- 0
# filter bad logger temp values
tracks$temp_mean[tracks$temp_mean > 20] <- NA

# Investigate dive histograms
plot.distributions(tracks, dive_cols)

# make distribtuion plots of steps, angles and dive param
data <- prepData.wrapper(tracks)
plot.distributions(data, c('step','angle','diveDuration_sum'))

##################
# HMM parameters #
##################
# Number of behavioural states
nbStates <- 2
# Probabiility distributions of data streams
dist = list(step="gamma", angle="vm", diveDuration_sum="gamma")
# state names
stateNames <- c('foraging', 'transit')
# Initial step distribution natural scale parameters
stepPar0 <- c(.1, .35, 2, 2)
# Initial angle distribution natural scale parameters
anglePar0 <- c(0, 0, 1, 20)
# Initial dive distribution natural scale parameters
divePar0 <- c(140, 20, 40, 20, 0, 0.6)
# Initial parmameter list
Par0 <- list(step=stepPar0, angle=anglePar0, diveDuration_sum=divePar0)
# Covariate formula
# covFormula <- ~ krig_indicator_low_10 + krig_indicator_medium_10 + krig_indicator_large_10 +
#   krig_indicator_extrem_10
covFormula <- ~ krig_salinity

##############
# Processing #
##############
# Produce a 2-state HMM model using momentuHMM
m <- fit.momentuHMM(tracks, nbStates, dist, Par0, covFormula, stateNames, 
                    standardiseCovs=T, estAngleMean=list(angle=TRUE))
beep()
saveRDS(m, model_out)
plot(m, plotStationary=T, plotCI=T, breaks=20, plotTracks=FALSE)
#par(ask=FALSE) # reset ask if quit early

########
# Maps #
########
plot.HMM.2state(m, type='combined_points', zoom=10, suffix='preyfield_paper')
plot.HMM.2state(m, type='split_points', zoom=10, suffix='preyfield_paper')
plot.HMM.2state(m, type='density', zoom=10, suffix='preyfield_paper')
# splitby survey
plot.HMM.2state.density.survey(m, tracks$krig_survey_id, zoom=10, suffix='preyfield_paper')

# plot HMM model foraging density over surveys
tracks$state <- viterbi(m)
# load the right surveys
surveys <- unique(tracks$krig_survey_id)
dbsAggs <- lapply(paste0('./kriging/output/agg/models/',
                         '3Dkrig_aggdb_Sv_mean_',surveys,
                         '.rds'), function(x) readRDS(x)@items)
track_split <- split(tracks, tracks$krig_survey_id)

# make the plots of top 30m
for (i in 1:length(track_split)){
  # tracks
  df <- lonlat2UTM.track(track_split[[i]])
  df$state <- as.factor(df$state)
  levels(df$state) <- stateNames
  df <- df[df$state == 'foraging',]
  # database
  db <- dbsAggs[[i]]
  # round top 10m
  db <- db[db$x3 <= 30 & db$Polygon,]
  db <- aggregate(db, by=list(db$x1, db$x2), FUN=mean)
  vmin <- min(db$Kriging.Sv_mean.estim, na.rm=T)
  vmax <- max(db$Kriging.Sv_mean.estim, na.rm=T)
  # make plot
  p <- ggplot(db, aes(x=x1, y=x2, fill=Kriging.Sv_mean.estim)) +
    geom_raster() +
    scale_fill_gradientn(colours=rev(cmocean('deep')(100)), 
                         limits=c(vmin, vmax), na.value="#A0A0A0") +
    geom_density_2d(aes(lon, lat, z=..level.., alpha=..level..), df, inherit.aes=F,
                    color='red') +
    ggtitle(surveys[i]) +
    scale_alpha(guide = 'none')
  print(p)
  invisible(readline(prompt=paste0("Showing ",surveys[i]," Press [enter] to continue")))
}

################
# Model checks #
################
# decode states and model data
tracks_out <- tracks
tracks_out$state <- as.factor(viterbi(m))
levels(tracks_out$state) <- stateNames
tracks_out[,c('step','angle')] <- m$data[,c('step','angle')]
# plot stream distributions
alpha=.4
bins=20
line_colors=c('#E69F00','#57B4E9')
p1 <- ggplot(tracks_out, aes(x=step, y=..density.., col=state)) +
  geom_histogram(fill=alpha('grey',alpha), col=alpha('darkgrey',alpha), bins=bins) +
  geom_density(size=.8) +
  scale_color_manual(values=line_colors) +
  theme(legend.position='none')
p2 <- ggplot(tracks_out, aes(x=angle, y=..density.., col=state)) +
  geom_histogram(fill=alpha('grey',alpha), col=alpha('darkgrey',alpha), bins=bins) +
  geom_density(size=.8) +
  scale_color_manual(values=line_colors)
p3 <- ggplot(tracks_out, aes(x=diveDuration_sum, y=..density.., col=state)) +
  geom_histogram(fill=alpha('grey',alpha), col=alpha('darkgrey',alpha), bins=bins) +
  geom_density(size=.8) +
  scale_color_manual(values=line_colors) +
  theme(legend.position='none')
# correaltion matrix of data streams
corr_matrix <- cor(tracks_out[,c('step','angle','diveDuration_sum')], use='pairwise.complete.obs')
p4 <- ggcorrplot(corr_matrix, colors=cmocean('balance')(5)[2:4], lab=T, type='upper')
# arrange
ggarrange(p1,p2,p3,p4)





