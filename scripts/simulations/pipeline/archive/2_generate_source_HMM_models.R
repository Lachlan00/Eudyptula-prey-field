# Archived (produce multiple models)
# New method produce a single model

# pilot control script for simulation scripts
# produce a set of simulated tracks
source('./scripts/eudyptula.R')
source('./scripts/HMM/HMM_functions.R')
source('./scripts/quality_control.R')
source("./scripts/visual.R")
# libraries
library(beepr)

# load data
tracks <- readRDS('./data/analysis_datasets/simulations/tracks_reg5m_simulationProcessed.rds')
# Reduce to Survey Periods
window.ls <- list(2,2,2,5,5)
names(window.ls) <- 2015:2019
tracks <- tracks.survey.time.filter(tracks, window=window.ls)
# clean gaps
tracks <- tracks.reg.filter.gaps(tracks)
# Final quality control checks
tracks <- quality.GPS.fix.count(tracks, min_threshold=70, filter=T, plotMap=F)
tracks <- quality.gap.ratio(tracks, max_threshold=1.2, filter=T, plotMap=F)

# quick look at steps and angles
plot.step.angle.dist(tracks)

####################
# Gemma Parameters #
####################
nbStates <- 2
dist = list(step="gamma", angle="vm")
stateNames <- c('foraging', 'transit')
stepPar0 <- c(.1, .35, 2, 2)
anglePar0 <- c(0, 0, 1, 40)
Par0 <- list(step=stepPar0, angle=anglePar0)
keepCols <- 'survey_id'

##############
# Model 2015 #
##############
mdat <- track.date.filter(tracks, '2015-01-01', '2016-01-01')
# Construct model
m.2015 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE), keepCols=keepCols)
beep()

##############
# Model 2016 #
##############
mdat <- track.date.filter(tracks, '2016-01-01', '2017-01-01')
# Construct model
m.2016 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE), keepCols=keepCols)
beep()

######################
# Lachlan Parameters #
######################
nbStates <- 2
dist = list(step="gamma", angle="vm")
stateNames <- c('foraging', 'transit')
stepPar0 <- c(.1, .5, 2, 2)
anglePar0 <- c(0, 0, 1, 12)
Par0 <- list(step=stepPar0, angle=anglePar0)

##############
# Model 2017 #
##############
mdat <- track.date.filter(tracks, '2017-01-01', '2018-01-01')
# Construct model
m.2017 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE), keepCols=keepCols)
beep()

##############
# Model 2018 #
##############
mdat <- track.date.filter(tracks, '2018-01-01', '2019-01-01')
# Construct model
m.2018 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE), keepCols=keepCols)
beep()

###########################
# 2019 (Gemma) Parameters #
###########################
nbStates <- 2
dist = list(step="gamma", angle="vm")
stateNames <- c('foraging', 'transit')
stepPar0 <- c(.1, .35, 2, 2)
anglePar0 <- c(0, 0, 1, 40)
Par0 <- list(step=stepPar0, angle=anglePar0)
keepCols <- 'survey_id'

##############
# Model 2019 #
##############
mdat <- track.date.filter(tracks, '2019-01-01', '2020-01-01')
# Construct model
m.2019 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE), keepCols=keepCols)
beep()

# all models done
beep('fanfare')

# Save models
saveRDS(m.2015, './data/simulations/HMM_models/HMM_2015.rds')
saveRDS(m.2016, './data/simulations/HMM_models/HMM_2016.rds')
saveRDS(m.2017, './data/simulations/HMM_models/HMM_2017.rds')
saveRDS(m.2018, './data/simulations/HMM_models/HMM_2018.rds')
saveRDS(m.2019, './data/simulations/HMM_models/HMM_2019.rds')

############
# PLOTTING #
############
plot(m.2019, breaks=20, plotTracks=T)
par(ask=FALSE)

plot.HMM.2state(m.2015, type='combined_points')
plot.HMM.2state(m.2015, type='density')

plot.HMM.2state(m.2016, type='combined_points')
plot.HMM.2state(m.2016, type='density')

plot.HMM.2state(m.2017, type='combined_points')
plot.HMM.2state(m.2017, type='density')

plot.HMM.2state(m.2018, type='combined_points')
plot.HMM.2state(m.2018, type='density')

plot.HMM.2state(m.2019, type='combined_points')
plot.HMM.2state(m.2019, type='density')

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

#####################################
# Attempts to unify model paramters #
#####################################

###TEST #######
# HMM parameters
nbStates <- 2
dist = list(step="gamma", angle="vm")
stateNames <- c('foraging', 'transit')
stepPar0 <- c(.1, .35, 2, 2)
anglePar0 <- c(0, 0, 1, 20)
Par0 <- list(step=stepPar0, angle=anglePar0)

##############
# Model 2016 #
##############
mdat <- track.date.filter(tracks, '2016-01-01', '2017-01-01')
# Construct model
m.2016 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE))
beep()
plot(m.2016, breaks=20, plotTracks=T)
# reset ask if quit early
par(ask=FALSE) 

##############
# Model 2019 #
##############
mdat <- track.date.filter(tracks, '2019-01-01', '2020-01-01')
# Construct model
m.2019 <- fit.momentuHMM(mdat, nbStates, dist, Par0, stateNames=stateNames,
                         estAngleMean=list(angle=TRUE))
beep()
plot(m.2019, breaks=20, plotTracks=T)
# reset ask if quit early
par(ask=FALSE) 

