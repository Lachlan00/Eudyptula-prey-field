# cumulative_prefield pilot tests
library(ggplot2)

source('./scripts/eudyptula.R')
source('./scripts/quality_control.R')

# load processed track data
tracks <- readRDS('./data/analysis_datasets/HMM/tracks_reg5m_processed.rds')
# clean gaps
tracks <- tracks.reg.filter.gaps(tracks)
# Final quality control checks
tracks <- quality.GPS.fix.count(tracks, min_threshold=50, filter=T)
# filter for surveys
tracks <- tracks.survey.time.filter(tracks, window=1)
# drop multiday tracks
tracks <- tracks.filter.multiday(tracks)
# drop 2019_S1 (only two tracks)
tracks <- tracks.filter.surveys(tracks, '2019_S1')

# Calculate cumulative prey
tracks <- tracks.cumsum.Sv(tracks)
# Add relative timestamps for each track
tracks$dtAEST_rel <- relatvie.timestamps(tracks$dtUTC)
# Filter out tracks that don't cover 8:00 - 15:00

# Make some pretty plots
ggplot(tracks, aes(x=dtAEST_rel, y=krig_Sv_mean_cumsum, col=id)) +
  geom_path(size=.3) +
  theme(legend.position = 'none') +
  facet_wrap(~survey_id)

# Break down inidviual surveys
survey <- '2017_S2'
ggplot(tracks.filter.surveys(tracks, survey, keep=T), 
       aes(x=dtAEST_rel, y=krig_Sv_mean_cumsum)) +
  geom_path(size=.4) +
  theme(legend.position = 'none') +
  geom_smooth(linetype='dashed', method='lm', size=.3) +
  facet_wrap(~id, scales = "free")
  
  
  
  
  