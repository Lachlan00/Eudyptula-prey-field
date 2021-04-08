# Simdata preperation

# prepare data for HMM modeling
library(beepr)
source("./scripts/HMM/HMM_functions.R")
source("./scripts/eudyptula.R")
source("./scripts/visual.R")
source("./scripts/quality_control.R")

# setup
input_fn <- './data/gps/penguin_tracks_all_REG-5m.rds'
output_fn <- './data/analysis_datasets/simulations/tracks_reg5m_simulationProcessed.rds'

# Processing
tracks <- readRDS(input_fn)
# filter tracks before 2015 (no kriging models)
tracks <- track.date.filter(tracks, t_start='2015-01-01')
# add survey id
tracks <- tracks.nearest.survey(tracks, survey_time_fn='data/transects/survey_times_simPaper.csv')
# only keep one survey id for each year (filter 2016_S1 and 2019_S1)
# tracks <- tracks.filter.surveys(tracks, c('2016_S1','2019_S1')) (Obsolete due to new survey time csv)
tracks <- rows.and.levels(tracks, idcol='survey_id')
tracks <- rows.and.levels(tracks)
# remove gaps by replacing constant speeds with NA
#tracks <- tracks.reg.filter.gaps(tracks)
# Using HMM model remove rafting periods near island
nbStates <- 2
stepPar0 <- c(0.1,0.5,0.1,0.5)
anglePar0 <- c(0,0,2,12)
tracks <- remove.rafting(tracks, nbStates, stepPar0, anglePar0, verbose = T)
# save output
saveRDS(tracks, output_fn)
beep()

