# Prepare cumulative prey analysis dataset
library(stringr)
library(zoo)
source('./scripts/simulations/simulate.R')
source('./scripts/eudyptula.R')
source('./scripts/quality_control.R')

###################################
# 1. Load real and simulated data #
###################################
# Load in real tracks and process
tracks.real <- readRDS('./data/analysis_datasets/simulations/tracks_reg5m_simulationProcessed.rds')
# Reduce to Survey Periods (if balanceTime = T window.ls is not used)
window.ls <- list(4,1,1,7,12)
names(window.ls) <- 2015:2019
tracks.real <- tracks.survey.time.filter(tracks.real, balanceTime = F, window=window.ls,
                                         survey_time_fn='./data/transects/survey_times_simPaper.csv')

tracksCount.report(tracks.real, tracks.real$survey_id)

# clean gaps
tracks.real <- tracks.reg.filter.gaps(tracks.real)
# Final quality control checks
# Drop low number of fixes (5th percentile)
# Drop high gap ratio  (95th percentile?)
min_threshold <- quantile(table(tracks.real$id[!is.na(tracks.real$lon)]), probs=.1)
max_threshold <- quantile(unlist(lapply(split(tracks.real, tracks.real$id),
                                        function(x) sum(is.na(x$lat))/sum(!is.na(x$lat)))), probs=.9)
tracks.real <- quality.GPS.fix.count(tracks.real, min_threshold=min_threshold, filter=T, plotMap=F)
tracks.real <- quality.gap.ratio(tracks.real, max_threshold=max_threshold, filter=T, plotMap=F)

tracks.real$speed <- NULL
tracks.real$type <- 'real'
tracks.real <- rows.and.levels(tracks.real)

tracksCount.report(tracks.real, tracks.real$survey_id)

# Load in simulations
tracks.sim <- readRDS('./data/simulations/sim_tracks.rds')
tracks.sim$id <- str_pad(tracks.sim$id, nchar(as.character(tracks.sim$id[nrow(tracks.sim)])),
                         pad = "0")
tracks.sim$id <- as.factor(tracks.sim$id)
tracks.sim$type <- 'sim'

# SAVE OUTPUT
saveRDS(tracks.real, './data/analysis_datasets/cumsum/tracks_real.rds')
saveRDS(tracks.sim, './data/analysis_datasets/cumsum/tracks_sim.rds')

# Below steps are now all done in the analysis script as we are now using a single 
# simulation set for evey survey so we need to recalcuate the cumulative Sv mean for
# each year

###########################################
# 2. Filter tracks that leave survey area #
###########################################
####(now done in analysis script)#####
# tracks <- inside.survey.zone(tracks.sim, threshold=10)
# # check filtering
# sim.report.tracksCount(tracks)

##########################################
# 3. Add environmental data to dataframe #
##########################################
####(now done in analysis script as needs to be repeated for each sim)#####
# Get covariate data for tracks
#tracks <- db.extract.covar.depth(tracks, depth=30, var='Sv_mean')

#########################################################################
# 4. Transform Sv_mean to linear values and interporlate missing values #
#########################################################################
####(now done in analysis script as needs to be repeated for each sim)#####
# # Transform Sv_mean
# tracks$Sv_linear <- Sv_mean.linear(tracks$krig_Sv_mean)
# # Interporlate missing values linearly
# tracks$Sv_linear <- na.approx(tracks$Sv_linear)

###############################
# 5. Calculate Cumlative Sums #
###############################
#tracks <- tracks.cumsum.Sv(tracks)



