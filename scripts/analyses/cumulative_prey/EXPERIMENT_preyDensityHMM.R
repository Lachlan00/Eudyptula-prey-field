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
window.ls <- list(2,2,2,5,10)
names(window.ls) <- 2015:2019
tracks <- tracks.survey.time.filter(tracks, window=window.ls)
# clean gaps
tracks <- tracks.reg.filter.gaps(tracks)
# Final quality control checks
tracks <- quality.GPS.fix.count(tracks, min_threshold=80, filter=T, plotMap=F)
tracks <- quality.gap.ratio(tracks, max_threshold=1.2, filter=T, plotMap=F)

# Filter qual number of tracks for each year
tracks <- equalise.tracks(tracks, tracks$survey_id, seed=200)
# report tracks by year (survey)
tracksCount.report(tracks, tracks$survey_id)
tracks.map(tracks, legend=F, zoom=10)

# quick look at steps and angles
plot.step.angle.dist(tracks)

# Add Sv mena to 30 m
tracks$krig_Sv_mean <- db.extract.covar.depth(tracks, 30, var='Sv_mean', method='mean')$krig_Sv_mean
tracks$krig_Sv_mean_linear <- Sv_mean.linear(tracks$krig_Sv_mean)

# add dive vars
tracks <- tracks.attach.diveVars()

# make covariate formula
covFormula <- ~ krig_Sv_mean_linear

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

######################
# Lachlan Parameters #
######################
nbStates <- 2
dist = list(step="gamma", angle="vm")
stateNames <- c('foraging', 'transit')
stepPar0 <- c(.1, .5, 2, 2)
anglePar0 <- c(0, 0, 1, 12)
Par0 <- list(step=stepPar0, angle=anglePar0)

###########################
# 2019 (Gemma) Parameters #
###########################
nbStates <- 2   
dist = list(step="gamma", angle="vm")
stateNames <- c('foraging', 'transit')
stepPar0 <- c(.1, .35, 2, 2)
anglePar0 <- c(0, 0, 1, 40)
divePar0 <- c(140, 20, 40, 20, 0, 0.6)
Par0 <- list(step=stepPar0, angle=anglePar0, diveDuration_sum=divePar0)
keepCols <- 'survey_id'

##############
# Model ALL #
##############
# Construct model
m.all <- fit.momentuHMM(tracks, nbStates, dist, Par0, stateNames=stateNames,
                        estAngleMean=list(angle=TRUE), keepCols=keepCols,
                        formula=covFormula)
beep()

# Save models
# saveRDS(m.all, './data/simulations/HMM_models/HMM.rds')
# saveRDS(tracks, './data/simulations/HMM_models/HMM_trackSource.rds')

############
# PLOTTING #
############
plot(m.all, breaks=20, plotTracks=T)
par(ask=FALSE)
plot.HMM.2state(m.all, type='combined_points')
plot.HMM.2state(m.all, type='density')

################################################
# SV mean transitional probability experiments #
################################################
