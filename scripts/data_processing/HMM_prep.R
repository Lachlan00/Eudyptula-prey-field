# prepare data for HMM modeling
library(beepr)
source("./scripts/HMM/HMM_functions.R")
source("./scripts/eudyptula.R")
source("./scripts/visual.R")
source("./scripts/quality_control.R")

# setup
input_fn <- './data/gps/penguin_tracks_all_REG-5m.rds'
output_fn <- './data/analysis_datasets/HMM/tracks_reg5m_processed.rds'

# Processing
tracks <- readRDS(input_fn)
# filter tracks before 2015 (no kriging models)
tracks <- track.date.filter(tracks, t_start='2015-01-01')
# add dive data
tracks <- tracks.attach.diveVars(tracks)
# get covariate data from kriging models
tracks <- tracks.attach.krig.covar(tracks, covars=c('temperature','salinity','Sv_mean','bathymetry'))
# add survey id
tracks <- tracks.nearest.survey(tracks)
# attach dive stats
dive_stats <- tracks_diveStats(tracks, return_dive_df=T)
tracks <- dive_stats[[1]]
dive_df <- dive_stats[[2]]
saveRDS(dive_df, './data/dive/dives_all.rds')
# remove gaps by replacing constant speeds with NA
#tracks <- tracks.reg.filter.gaps(tracks)
# Usinbg HMM model remove rafting periods near island
nbStates <- 2
stepPar0 <- c(0.1,0.5,0.1,0.5)
anglePar0 <- c(0,0,2,12)
tracks <- remove.rafting(tracks, nbStates, stepPar0, anglePar0)

# save output
saveRDS(tracks, output_fn)
beep()

