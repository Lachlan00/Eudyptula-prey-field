# Tables
source('scripts/eudyptula.R')

#################
# Penguin table #
#################
# Load penguins
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
# Filter tracks that leave the survey area like in previous analysis
tracks.real <- inside.survey.zone(tracks.real, threshold=10, 
                                  plot.map=F, plot.title='Real Tracks')

# Tracks per year
tracks.real$year <- year(tracks.real$dtUTC)
tracksCount.report(tracks.real, tracks.real$year)
# Individuals per year
tracks.ind <- tracks.real[!duplicated(str_sub(tracks.real$id,1,5)),]
table(tracks.ind$year)
# No fixes
table(tracks.real$year)
# Raw fixes
tracks.raw <- readRDS('./data/gps/penguin_tracks_all_CLEAN.rds')
tracks.raw <- tracks.raw[tracks.raw$id %in% as.character(tracks.ind$id),]
table(year(tracks.raw$dtUTC))

# mean time between fixes
tracks.raw <- tracks.time(tracks.raw)
aggregate(tracks.raw$time, by=list(year(tracks.raw$dtUTC)), function(x)  median(x, na.rm = T))

# Deployment window
aggregate(tracks.real$dtUTC, by=list(tracks.real$year), function(x)  min(x, na.rm = T))
aggregate(tracks.real$dtUTC, by=list(tracks.real$year), function(x)  max(x, na.rm = T))

# load penguin dive data
dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
# dives longer than 5 seconds
longCutoff <- 5
dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > longCutoff,]
dive.data.long$year <- year(dive.data.long$start)
dive.data.long <- dive.data.long[dive.data.long$id %in% as.character(tracks.ind$id),]
table(dive.data.long$year)

# Sex Distribution
tracks.ind <- sex.factor(tracks.ind)
aggregate(tracks.ind$sex, by=list(tracks.ind$year), table)

# Stats
tracks.real <- rows.and.levels(tracks.real)
sum.stats <- peng.summary.stats(tracks.real)
aggregate(sum.stats[,c('duration', 'displacement_max', 'distance_total')],
          by=list(sum.stats$year), mean)
aggregate(sum.stats[,c('duration', 'displacement_max', 'distance_total')],
          by=list(sum.stats$year), sd)

