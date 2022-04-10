library(ks)
library(ggplot2)
source('scripts/eudyptula.R')
source('scripts/kde_percent.R')

# load acoustics
surveys = c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
krig.dfs <- load.krig.dbs('Sv_mean', surveys, return_items = T)
for (i in 1:length(krig.dfs)){
  krig.dfs[[i]]$year <- (2015:2019)[i]
}
krig.dfs <- do.call(rbind, krig.dfs)
krig.dfs <- krig.dfs[krig.dfs$x3 <= 30,]
krig.dfs$linear <- Sv_mean.linear(krig.dfs$Kriging.Sv_mean.estim)

# Calculate Sv mean in core range
mod.HMM <- readRDS('./data/simulations/HMM_models/HMM.rds')
m.tracks <- mod.HMM$data
m.tracks$state <- viterbi(mod.HMM)
names(m.tracks)[c(1,4,5,6,7)] <- c('id','lon','lat','survey','state')
# Filter tracks that leave the survey area like in previous analysis
m.tracks <- inside.survey.zone(m.tracks, threshold=20, use.percentage = T,
                               plot.map=F, plot.title='Real Tracks')
# get only foraging state
m.tracks <- m.tracks[m.tracks$state == 1,]
m.tracks$survey <- as.factor(substr(m.tracks$survey,1,4))
kde.df <- KDE_percentage(m.tracks[complete.cases(m.tracks),c('lon','lat')],
                         levels=c(kde.per), facet=tracks[complete.cases(m.tracks),c('survey')])
kde.df$survey <- kde.df$facet
# find what levels are for each survey
kernel_data.ls <- split(kde.df, kde.df$survey)
# For each kernel data set extract the outside border and create
# shape files from them
kernel_poly <- list()
for (i in 1:length(kernel_data.ls)){
  kdata <- kernel_data.ls[[i]]
  # Split and make into shape files
  kdata.ls <- split(kdata, kdata$group)
  # make polygons
  for (j in 1:length(kdata.ls)){
    kdata.ls[[j]] <- Polygon(cbind(kdata.ls[[j]]$x, kdata.ls[[j]]$y))
  }
  kernel_poly[[i]] <- SpatialPolygons(list(Polygons(kdata.ls, 1)))
}
names(kernel_poly) <- 2015:2019

# Extract Sv from polygons
names(krig.dfs)[2:3] <- c('lon', 'lat')
krig.dfs <- UTM2lonlat.df(krig.dfs)

# calc over
krig.dfs.ls <- split(krig.dfs, krig.dfs$year)
for (i in 1:length(krig.dfs.ls)){
  dat <- data.frame(lon = krig.dfs.ls[[i]]$lon,
                    lat = krig.dfs.ls[[i]]$lat)
  coordinates(dat) <- ~ lon + lat
  krig.dfs.ls[[i]]$over <- !is.na(over(dat, kernel_poly[[i]]))
}
krig.dfs <- do.call(rbind, krig.dfs.ls)

bool <- krig.dfs$over
krig.dfs$over[bool] <- 'Inside'
krig.dfs$over[!bool] <- 'Outside'

# calc grand means
grand.mean <- mean(krig.dfs$linear)
krig.dfs$diff.linear <- krig.dfs$linear - grand.mean

# Now do some kind of analysis
ggplot() +
  geom_boxplot(data=krig.dfs, outlier.shape = NA,
               mapping=aes(x=factor(year),
                           y=diffs))+#Kriging.Sv_mean.estim,
                           #fill=factor(over))) +
               #lims(y=c(-70,-40)) +
  lims(y=c(-3e-6, 1e-5)) +
  ylab("Sv mean") +
  xlab(element_blank()) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(text = element_text(size=18)) +
  labs(fill='75% UD') +
  geom_hline(yintercept=0, linetype='dashed')


ggplot(krig.dfs, aes(x=factor(year), y=diff.linear)) +
  geom_hline(yintercept = 0, color='darkred', linetype='dashed', size=.7) +
  geom_boxplot(fill='grey', alpha=.6, outlier.shape = NA) +
  labs(x=NULL, y='Sv mean anomaly (linear)', fill="Sex") +
  theme_bw() + grids(linetype = "dashed") +
  theme(text = element_text(size=18)) +
  theme(legend.position = "none") +
  #lims(y=c(-3e-6, 1e-5)) +
  geom_label(data=text.Sv, mapping=aes(x=x, y=y, label=label),
             inherit.aes = F, fill='grey', alpha=.6)


