# Data checks to make sure 2017 differences are not a sampling error
source('./scripts/eudyptula.R')
# Read in data in the same manner as the analysis

# Analysis paramters
depth = 30 # How deep to get the mean of Sv mean
weight.depth = T # Weight the depth aggregation by penguin diving behaviour
boundary_thresh = 10 # Filter tracks with points > N outside survey boundary
track_regTime = 5 # Time between fixes in minutes (match to tracks)
bootstrap_N = 100000 # Iterations for bootstrap test
equalise_tracks = F # Reduce real tracks so equal number for each year?
set.seed(420)

# Load the data
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
tracks.sim <- readRDS('./data/analysis_datasets/cumsum/tracks_sim.rds')

# Filter tracks that leave the survey area
tracks.real <- inside.survey.zone(tracks.real, threshold=boundary_thresh, 
                                  plot.map=F, plot.title='Real Tracks')
tracks.sim <- inside.survey.zone(tracks.sim, threshold=boundary_thresh,
                                 plot.map=F, plot.title='Simulated Tracks')
tracksCount.report(tracks.real, tracks.real$survey_id)

# Equalise tracks
if (equalise_tracks){
  tracks.real <- equalise.tracks(tracks.real, tracks.real$survey_id)
}

# Load the survey times
survey.times <- read.csv('./data/transects/survey_times_simPaper.csv')
survey.times$start_UTC <- as.POSIXct(survey.times$start_UTC)
survey.times$end_UTC <- as.POSIXct(survey.times$end_UTC)
survey.times <- data.frame(datetime=c(survey.times$start_UTC, survey.times$end_UTC))

# Get mean track times
track.times <- aggregate(tracks.real$dtUTC, by=list(tracks.real$id), mean)
names(track.times)  <- c('id', 'datetime')
# Plot the stuff
ggplot() +
  geom_point(mapping=aes(x=datetime, y=1), data=track.times) +
  geom_line(mapping=aes(x=datetime, y=2), data=survey.times, col='red') +
  facet_wrap(~year(datetime), scale='free_x', ncol=1) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
