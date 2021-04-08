# Dive data
source('scripts/eudyptula.R')
source('scripts/HMM/HMM_functions.R')
source('scripts/visual.R')

dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')[[1]]
# Tracks analysed
tracks <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')

# Filter dives
dive.data  <- dive.data[as.character(dive.data$id) %in% as.character(unique(tracks$id)),]
# Keep only ones with location
dive.data <- dive.data[complete.cases(dive.data[,]),]
# Add Sv_mean
dive.data <- tracks.attach.krig.covar(dive.data, covars='Sv_mean', 
                                      depth_colname = 'diveDepth_max')

# Quick plots
ggplot(dive.data, aes(x=krig_Sv_mean, y=diveDuration_sum)) +
  geom_point(alpha=.1) +
  geom_smooth(method='loess') +
  facet_wrap(~krig_survey_id)
  
# Map
names(dive.data)[5] <- 'survey_id'
tracks.map(dive.data[complete.cases(dive.data),], colfactor='krig_Sv_mean', 
           facet_survey=T, plot_line=F, zoom=10)
