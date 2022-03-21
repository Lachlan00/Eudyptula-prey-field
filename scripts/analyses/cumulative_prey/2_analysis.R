library(ggplot2)
library(zoo)
library(lazyWeave)
library(beepr)
source('./scripts/eudyptula.R')

# Analysis paramters
depth = 30 # How deep to get the mean of Sv mean
weight.depth = T # Weight the depth aggregation by penguin diving behaviour
boundary_thresh = 100 # Filter tracks with points > N outside survey boundary
thresh.use.percentage = T # Above filter becomes percentage of points outside survey area
track_regTime = 5 # Time between fixes in minutes (match to tracks)
bootstrap_N = 100000 # Iterations for bootstrap test
equalise_tracks = F # Reduce real tracks so equal number for each year?

set.seed(420)

####################
# 1. Preprocessing #
####################
# Load the data
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
tracks.sim <- readRDS('./data/analysis_datasets/cumsum/tracks_sim.rds')

# Filter tracks that leave the survey area
tracks.real <- inside.survey.zone(tracks.real, threshold=boundary_thresh, add.inZone.col=T, 
                                  plot.map=F, plot.title='Real Tracks', use.percentage=thresh.use.percentage)
tracks.sim <- inside.survey.zone(tracks.sim, threshold=boundary_thresh, add.inZone.col=T,
                                 plot.map=F, plot.title='Simulated Tracks', use.percentage=thresh.use.percentage)
# make a secondary set of simulations that don't use the area to the north
# that wasn't sampled in 2015
# filter all simulations that went north of -36.235 lat
tracks.sim.2015 <- split(tracks.sim, factor(tracks.sim$id))
for (i in 1:length(tracks.sim.2015)){
  if (max(tracks.sim.2015[[i]]$lat) > -36.235){
    tracks.sim.2015[[i]] <- NA
  }
}
tracks.sim.2015 <- tracks.sim.2015[!is.na(tracks.sim.2015)]
tracks.sim.2015 <- do.call(rbind, tracks.sim.2015)

# Report track info
tracksCount.report(tracks.real, tracks.real$survey_id)

# Equalise tracks
if (equalise_tracks){
  tracks.real <- equalise.tracks(tracks.real, tracks.real$survey_id)
}

# If using a weighted mean load dive data
if (weight.depth){
  # load penguin dive data
  dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
  # dives longer than 5 seconds
  longCutoff <- 5
  dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > longCutoff,]
  dive.data.long$year <- year(dive.data.long$start)
  dive.data.long <- dive.data.long[as.character(dive.data.long$id) %in% 
                                     as.character(tracks.real$id),]
  dive.data.long <- tracks.nearest.survey(dive.data.long, dtUTCcol='start')
  dive.data.long <- rows.and.levels(dive.data.long, idcol='survey_id')
  # Split and calculate the distribution
  dive.dist <- data.frame(survey_id=dive.data.long$survey_id,
                          depth=as.integer(round(dive.data.long$depth_max)))
  dive.dist.ls <- split(dive.dist, dive.dist$survey_id)
  for (i in 1:length(dive.dist.ls)){
    res <- data.frame(depth=0:max(dive.dist$depth),
                      freq=0)
    tab <- table(dive.dist.ls[[i]]$depth)
    res[as.integer(names(tab))+1,'freq'] <- tab
    res$survey_id <- dive.dist.ls[[i]]$survey_id[1]
    dive.dist.ls[[i]] <- res
  }
  dive.weights <- do.call(rbind, dive.dist.ls)
  dive.weights <- rows.and.levels(dive.weights, idcol='survey_id')
}

# create an index column
tracks.real <- tracks.add.indexes(tracks.real)
tracks.sim <- tracks.add.indexes(tracks.sim)
tracks.sim.2015 <- tracks.add.indexes(tracks.sim.2015)

# Get covariate data for  real tracks
if (weight.depth){
  tracks.real <- db.extract.covar.depth(tracks.real, method='weight', var='Sv_mean', 
                                        weight.df=dive.weights)
} else {
  tracks.real <- db.extract.covar.depth(tracks.real, depth=depth, 
                                        var='Sv_mean', weight.df=dive.weights)
}
# Convert to linear
tracks.real$Sv_linear <- Sv_mean.linear(tracks.real$krig_Sv_mean)
# Interporlate missing values linearly
# Update! Don't interporalte, just remove
#tracks.real$Sv_linear <- na.approx(tracks.real$Sv_linear)
tracks.real <- tracks.real[!is.na(tracks.real$Sv_linear),]
# Remove points outside of zone
tracks.real <- tracks.real[tracks.real$inZone,]

# Calculate cumulatice preyfield
tracks.real <- tracks.cumsum.Sv(tracks.real)
track.real.counts <- tracksCount.report(tracks.real, split.var=tracks.real$survey_id)

# Split tracks real by year
tracks.real.ls <- split(tracks.real, tracks.real$survey_id)

############
# Analysis #
############
# list for plots
plt.lines <- list()
plt.box <- list()
plt.boot.hist <- list()
stat.tests <- list()
analysis.ls <- list()
# Loop through surveys of tracks real and compare agianst simulations
for (i in 1:length(tracks.real.ls)){
  # Get the real tracks for the survey
  survey <- names(tracks.real.ls)[i]
  message('\nAnalysing ',survey,'..')
  tracks.real <- tracks.real.ls[[i]]
  # Process the simulations
  # Get Sv mean (convert to UTM only on first pass)
  # For other sims that is i == 2
  if (weight.depth) {
    if (survey == '2015_S1'){
      tracks.sim.2015 <- db.extract.covar.depth(tracks.sim.2015, method='weight', 
                                           convert2UTM=T,  
                                           var='Sv_mean', survey=survey,
                                           weight.df = dive.weights)
    } else {
      tracks.sim <- db.extract.covar.depth(tracks.sim, method='weight', 
                                           convert2UTM=ifelse(i==2,T,F),  
                                           var='Sv_mean', survey=survey,
                                           weight.df = dive.weights)
    }
  } else {
    if (survey == '2015_S1'){
      tracks.sim.2015 <- db.extract.covar.depth(tracks.sim.2015, depth=depth, 
                                           convert2UTM=T,  
                                           var='Sv_mean', survey=survey)
    } else {
      tracks.sim <- db.extract.covar.depth(tracks.sim, depth=depth, 
                                           convert2UTM=ifelse(i==2,T,F),  
                                           var='Sv_mean', survey=survey)
    }
  }
  # set what sims we are using
  if (survey == '2015_S1'){
    track.sims.active <- tracks.sim.2015
  } else {
    track.sims.active <- tracks.sim
  }
  # Convert to linear
  track.sims.active$Sv_linear <- Sv_mean.linear(track.sims.active$krig_Sv_mean)
  # Remove points outside of zone
  track.sims.active <- track.sims.active[track.sims.active$inZone,]
  # Calculate cumulatice preyfield
  track.sims.active <- tracks.cumsum.Sv(track.sims.active)
  
  # Now create an analysis dataframe
  df.analysis <- rbind(tracks.real[,c('id','rank','type','Sv_linear_cumsum')],
                       track.sims.active[,c('id','rank','type','Sv_linear_cumsum')])
  df.analysis <- rows.and.levels(df.analysis)
  
  # Now for each track, calulcate the consumption rate per hour
  df.consumption <- split(df.analysis, df.analysis$id)
  type.ls <- list()
  for (j in 1:length(df.consumption)){
    type.ls[[j]] <- df.consumption[[j]]$type[1]
    df.consumption[[j]] <- df.consumption[[j]]$Sv_linear_cumsum[
      nrow(df.consumption[[j]])]/(nrow(df.consumption[[j]])*track_regTime/60)
  }
  df.consumption <- data.frame(id=names(df.consumption),
                               type=unlist(type.ls),
                               consumption=unlist(df.consumption))
  df.consumption <- rows.and.levels(df.consumption)
  # add to analysis set
  analysis.ls[[i]] <- list(cumulative=df.analysis, consumption=df.consumption)
  
  # Perform the anlysis
  # __Visual analysis__ #
  # Line plot
  plt.lines[[i]] <- ggplot(analysis.ls[[i]]$cumulative, 
                           aes(x=rank, y=Sv_linear_cumsum, col=type)) +
    geom_smooth(method='lm') +
    ggtitle(survey) +
    labs(y='Cumulative prey encounters', x=NULL) +
    theme(legend.position = 'none')
  # Boxplot
  plt.box[[i]] <- ggplot(analysis.ls[[i]]$consumption,
                         aes(x=type, y=consumption, fill=type)) +
    geom_boxplot() +
    ggtitle(survey) +
    labs(y='Prey encounter per hour', x=NULL) +
    theme(legend.position = 'none')

  # Stats
  # Standard lm model
  stat.tests[[i]] <- list(lm=lm(formula = consumption ~ type, analysis.ls[[i]]$consumption))
  # Bootstrapped test
  # Sample sims N times and see how often then mean 
  # consumption rate is lower than the real tracks
  bootstrapped.sim.means <- rep(0, bootstrap_N)
  for (j in 1:bootstrap_N) {
    bootstrapped.sim.means[j] <- mean(sample(analysis.ls[[i]]$consumption$consumption[
                                                analysis.ls[[i]]$consumption$type == 'sim'],
                                             size=track.real.counts$count[i]))
  }
  analysis.ls[[i]]$bootstrapped.sim.means <- bootstrapped.sim.means
  
  # calculate proportion of real > sim
  stat.tests[[i]]$bootstrap.prop <- sum(bootstrapped.sim.means < mean(analysis.ls[[i]]$consumption$consumption[
                                      analysis.ls[[i]]$consumption$type == 'real'])) / bootstrap_N
  # make a histogram
  plt.boot.hist[[i]] <- ggplot() + 
    geom_histogram(mapping=aes(x=bootstrapped.sim.means), bins=20, 
                   col='grey', fill='darkgrey') +
    geom_vline(xintercept=mean(analysis.ls[[i]]$consumption$consumption[
      analysis.ls[[i]]$consumption$type == 'real']), col='red', linetype='solid', size=1) +
    geom_vline(xinterScept=mean(analysis.ls[[i]]$consumption$consumption[
      analysis.ls[[i]]$consumption$type == 'real']) - 
        sd(analysis.ls[[i]]$consumption$consumption[
        analysis.ls[[i]]$consumption$type == 'real']), col='red', linetype='dashed', size=.5) +
    geom_vline(xintercept=mean(analysis.ls[[i]]$consumption$consumption[
      analysis.ls[[i]]$consumption$type == 'real']) + 
        sd(analysis.ls[[i]]$consumption$consumption[
          analysis.ls[[i]]$consumption$type == 'real']), col='red', linetype='dashed', size=.5) +
    labs(x='Simulated means of prey encounter rate', y='Count') +
    ggtitle(paste0(survey,' (',round(stat.tests[[i]]$bootstrap.prop,2),')'))
}

# Present plots
# ggarrange(plotlist=plt.lines)
# ggarrange(plotlist=plt.box)
print(ggarrange(plotlist=plt.boot.hist))

# Present stats
# lm
print(data.frame(survey=names(tracks.real.ls),
                 pval=pvalString(unlist(lapply(stat.tests, 
                                    function(m) summary(m$lm)$coefficients[2,4])))))
# bootstrap
print(data.frame(survey=names(tracks.real.ls),
                 prop=unlist(lapply(stat.tests, 
                                    function(m) m$bootstrap.prop))))

beep()
# Save data for better plotting later
# saveRDS(analysis.ls, './data/analysis_datasets/cumsum/simPaper_analysis.rds')
# saveRDS(stat.tests, './data/analysis_datasets/cumsum/simPaper_stats.rds')
