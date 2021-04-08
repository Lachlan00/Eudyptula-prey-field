# Post hoc analysis to see why 2017 is different
library(ggplot2)
library(ggcorrplot)
library(cmocean)
library(anchors)
library(mgcv)
library(visreg)
source('scripts/eudyptula.R')
source('scripts/visual.R')
source('scripts/IMOS.R')
source('scripts/HMM/HMM_functions.R')

# Paramaters
depth = 30 # Depth to extract Sv_mean data

# load databases
surveys <- c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
dbs.Svmn <- load.krig.dbs(var='Sv_mean', surveys=surveys, return_items=T)
dbs.temp <- load.krig.dbs(var='temperature', surveys=surveys, return_items=T)
dbs.salt <- load.krig.dbs(var='salinity', surveys=surveys, return_items=T)

# Extract the top water of water column to depth
dbs.Svmn <- lapply(dbs.Svmn, function(db) db[db$x3 <= depth & db$Polygon,])
dbs.temp <- lapply(dbs.temp, function(db) db[db$x3 <= depth & db$Polygon,])
dbs.salt <- lapply(dbs.salt, function(db) db[db$x3 <= depth & db$Polygon,])

# Take the mean over x3 (depth)
# Comment out these lines to avoid doing the mean for water column
dbs.Svmn.agg <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x1, db$x2), FUN=mean))
dbs.temp.agg <- lapply(dbs.temp, function(db) aggregate(db, by=list(db$x1, db$x2), FUN=mean))
dbs.salt.agg <- lapply(dbs.salt, function(db) aggregate(db, by=list(db$x1, db$x2), FUN=mean))

# Aggregate in 2 dimesnions
dbs.Svmn.vert <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=mean))
dbs.temp.vert <- lapply(dbs.temp, function(db) aggregate(db, by=list(db$x3), FUN=mean))
dbs.salt.vert <- lapply(dbs.salt, function(db) aggregate(db, by=list(db$x3), FUN=mean))

# Add survey name to each database set
for (i in 1:length(surveys)){
  dbs.Svmn[[i]]$survey <- surveys[i]
  dbs.temp[[i]]$survey <- surveys[i]
  dbs.salt[[i]]$survey <- surveys[i]
  dbs.Svmn.agg[[i]]$survey <- surveys[i]
  dbs.temp.agg[[i]]$survey <- surveys[i]
  dbs.salt.agg[[i]]$survey <- surveys[i]
  dbs.Svmn.vert[[i]]$survey <- surveys[i]
  dbs.temp.vert[[i]]$survey <- surveys[i]
  dbs.salt.vert[[i]]$survey <- surveys[i]
}

#########################
# Vertical distribution #
#########################
# Add if survey 2017_S1
for (i in 1:length(surveys)){
  dbs.Svmn.vert[[i]]$survey2017 <- dbs.Svmn.vert[[i]]$survey[1] == '2017_S2'
  dbs.temp.vert[[i]]$survey2017 <- dbs.temp.vert[[i]]$survey[1] == '2017_S2'
  dbs.salt.vert[[i]]$survey2017 <- dbs.salt.vert[[i]]$survey[1] == '2017_S2'
}

# __With all surveys colored__ #
# SV_mean
pSvmn <- ggplot(do.call(rbind, dbs.Svmn.vert), 
       aes(x=x3, y=Kriging.Sv_mean.estim, col=survey)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  coord_flip() +
  ggtitle('Sv mean')

# Temperature
ptemp <- ggplot(do.call(rbind, dbs.temp.vert), 
       aes(x=x3, y=Kriging.temperature.estim, col=survey)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  coord_flip() +
  ggtitle('Temperature')

# Salinity
psalt <- ggplot(do.call(rbind, dbs.salt.vert), 
       aes(x=x3, y=Kriging.salinity.estim, col=survey)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  coord_flip() +
  ggtitle('Salinity')

ggarrange(plotlist = list(pSvmn, ptemp, psalt), common.legend=T, legend='right')
  
# __With 2017 colored__ #
# SV_mean
pSvmn <- ggplot(do.call(rbind, dbs.Svmn.vert), 
                aes(x=x3, y=Kriging.Sv_mean.estim, group=survey, col=survey2017)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  ggtitle('Sv mean')

# Temperature
ptemp <- ggplot(do.call(rbind, dbs.temp.vert), 
                aes(x=x3, y=Kriging.temperature.estim, group=survey, col=survey2017)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  ggtitle('Temperature')

# Salinity
psalt <- ggplot(do.call(rbind, dbs.salt.vert), 
                aes(x=x3, y=Kriging.salinity.estim, group=survey, col=survey2017)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  scale_color_manual(values=c('grey','#F8766D')) +
  coord_flip() +
  ggtitle('Salinity')

ggarrange(plotlist = list(pSvmn, ptemp, psalt), common.legend=T, legend='right')

##########################################
# Look at correlations between variables #
##########################################
# Make a master dataframe
dbdf <- do.call(rbind, dbs.Svmn.agg)[,c('survey','Kriging.Sv_mean.estim')]
names(dbdf)[2] <- 'Svmn'
row.names(dbdf) <- NULL
dbdf$temp <- do.call(rbind, dbs.temp.agg)[,'Kriging.temperature.estim']
dbdf$salt <- do.call(rbind, dbs.salt.agg)[,'Kriging.salinity.estim']
survey_rows <- dbdf$survey
dbdf$survey <- NULL
cor.ls <- split(dbdf, survey_rows)
# Correlate
cor.ls <- lapply(cor.ls, function(df) cor(df))

# Plot
p.cor.ls <- list()
for (i in 1:length(surveys)){
  p.cor.ls[[i]] <- ggcorrplot(cor.ls[[i]], 
                              colors=c('#0b3487','#ffffff','#87130b'), 
                              lab=T, type='upper') +
    ggtitle(surveys[i])
}
ggarrange(plotlist=p.cor.ls, common.legend=T, legend='right')

########
# Maps #
########
# __SV mean__ #
p.ls <- list()
vmin <- min(unlist(lapply(dbs.Svmn.agg, function(df) min(df$Kriging.Sv_mean.estim))))
vmax <- max(unlist(lapply(dbs.Svmn.agg, function(df) max(df$Kriging.Sv_mean.estim))))
for (i in 1:length(surveys)){
  p.ls[[i]] <- ggplot(data=dbs.Svmn.agg[[i]],
         aes(x=x1/1000, y=x2/1000,fill=Kriging.Sv_mean.estim)) +
    geom_raster() +
    scale_fill_gradientn(colours=rev(cmocean('deep')(100)),na.value="#A0A0A0",
                         limits=c(vmin, vmax)) +
    labs(fill='Sv mean') +
    ggtitle(paste(surveys[i])) +
    xlab("UTM Easting (km)") +
    ylab("UTM Northing (km)") +
    theme(legend.key.height = unit(2, "cm"))
}
ggarrange(plotlist=p.ls, common.legend=T, legend='right')

# __Temperature__ #
p.ls <- list()
vmin <- min(unlist(lapply(dbs.temp.agg, function(df) min(df$Kriging.temperature.estim))))
vmax <- max(unlist(lapply(dbs.temp.agg, function(df) max(df$Kriging.temperature.estim))))
for (i in 1:length(surveys)){
  p.ls[[i]] <- ggplot(data=dbs.temp.agg[[i]],
                      aes(x=x1/1000, y=x2/1000,fill=Kriging.temperature.estim)) +
    geom_raster() +
    scale_fill_gradientn(colours=cmocean('thermal')(100), na.value="#A0A0A0",
                         limits=c(vmin, vmax)) +
    labs(fill='Temperature (°C)') +
    ggtitle(paste(surveys[i])) +
    xlab("UTM Easting (km)") +
    ylab("UTM Northing (km)") +
    theme(legend.key.height = unit(2, "cm"))
}
ggarrange(plotlist=p.ls, common.legend=T, legend='right')

# __Salinity__ #
p.ls <- list()
vmin <- min(unlist(lapply(dbs.salt.agg, function(df) min(df$Kriging.salinity.estim))))
vmax <- max(unlist(lapply(dbs.salt.agg, function(df) max(df$Kriging.salinity.estim))))
for (i in 1:length(surveys)){
  p.ls[[i]] <- ggplot(data=dbs.salt.agg[[i]],
                      aes(x=x1/1000, y=x2/1000,fill=Kriging.salinity.estim)) +
    geom_raster() +
    scale_fill_gradientn(colours=cmocean('haline')(100), na.value="#A0A0A0",
                         limits=c(vmin, vmax)) +
    labs(fill='Salinity (PSU)') +
    ggtitle(paste(surveys[i])) +
    xlab("UTM Easting (km)") +
    ylab("UTM Northing (km)") +
    theme(legend.key.height = unit(2, "cm"))
}
ggarrange(plotlist=p.ls, common.legend=T, legend='right')

##############################
# Tracks against simulations #
##############################
# Load the data
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
tracks.sim <- readRDS('./data/analysis_datasets/cumsum/tracks_sim.rds')

# Filter tracks that leave the survey area like in previous analysis
tracks.real <- inside.survey.zone(tracks.real, threshold=10, 
                                  plot.map=F, plot.title='Real Tracks')
tracks.sim <- inside.survey.zone(tracks.sim, threshold=10,
                                 plot.map=F, plot.title='Simulated Tracks')

# __Make a map of real tracks__ #
tracks.map(tracks.real, colfactor='survey_id', zoom=10, 
           plot_line=F, facet_survey=T, legend=F)

# __Make a map of 2017 and sim tracks__ #
data <- tracks.real[tracks.real$survey_id == '2017_S2', c('lon','lat','type')]
data <- rbind(tracks.sim[,c('lon','lat','type')], data)
tracks.map(data,
           colfactor='type', zoom=10, 
           plot_line=F, legend=F)

################################
# __Real tracks over Sv mean__ #
################################
# Convert real tracks to UTM
tracks.real.UTM <- na.sub.FUN.df(tracks.real, c('lon','lat'), lonlat2UTM.df)
# Map databases
p.ls <- list()
for (i in 1:length(surveys)){
  p.ls[[i]] <- ggplot(data=dbs.Svmn.agg[[surveys[i]]],
              aes(x=x1/1000, y=x2/1000,fill=Kriging.Sv_mean.estim)) +
    geom_raster() +
    scale_fill_gradientn(colours=rev(cmocean('deep')(100)), na.value="#A0A0A0") +
    labs(fill='Sv mean') +
    xlab("UTM Easting (km)") +
    ylab("UTM Northing (km)") +
    theme(legend.key.height = unit(2, "cm"))
  # Add tracks
  p.ls[[i]] <- p.ls[[i]] + 
    geom_point(mapping=aes(x=lon/1000, y=lat/1000), 
                 data=tracks.real.UTM[tracks.real.UTM$survey_id == surveys[i],], 
               inherit.aes=F, alpha=0.08, col='red')
}
ggarrange(plotlist=p.ls, common.legend=T, legend='right')

##########################
# Downlaod survey images #
##########################
download = FALSE # Set false if already completed
if (download){
  survey.times <- load.survey.times(process=F)
  for (i in 1:length(surveys)){
    daterange <- do.call(c, as.list(survey.times[survey.times$id == surveys[i], c(3,4)]))
    IMOS.OceanCurrent.download('SST', 'SNSW', daterange, prefix=surveys[i])
    IMOS.OceanCurrent.download('chla', 'SNSW', daterange, prefix=surveys[i])
    IMOS.OceanCurrent.download('SST', 'Syd-Hob', daterange, prefix=surveys[i])
    IMOS.OceanCurrent.download('chla', 'Syd-Hob', daterange, prefix=surveys[i])
  }
}

##########################
# Sv mean center of mass #
##########################
# Calculate the center of mass for kriging maps
p.ls <- list()
cgi.ls <- list()
for (i in 1:length(surveys)){
  # Calc cgi (Sv mean must be linear)
  cgi.ls[[i]] <- calc.CGI(dbs.Svmn.agg[[i]]$x1, dbs.Svmn.agg[[i]]$x2, 
                          Sv_mean.linear(dbs.Svmn.agg[[i]]$Kriging.Sv_mean.estim))
  # Map center of gravity
  p.ls[[i]] <- ggplot(data=dbs.Svmn.agg[[i]],
         aes(x=x1/1000, y=x2/1000,fill=Kriging.Sv_mean.estim)) +
    geom_raster() +
    scale_fill_gradientn(colours=rev(cmocean('deep')(100)),na.value="#A0A0A0") +
    labs(fill='Sv mean') +
    ggtitle(paste(surveys[i])) +
    xlab("UTM Easting (km)") +
    ylab("UTM Northing (km)") +
    theme(legend.key.height = unit(2, "cm"))
  # Add CGI
  p.ls[[i]] <- p.ls[[i]] +
    geom_line(data=cgi.ls[[i]][cgi.ls[[i]]$name != 'centre_gravity',],
              aes(x=x/1000,y=y/1000, group=name), inherit.aes=F, 
              col='red', linetype='dashed', size=.35) +
    geom_point(data=cgi.ls[[i]][cgi.ls[[i]]$name == 'centre_gravity',],
               mapping=aes(x=x/1000, y=y/1000), inherit.aes=F,
               col='red', size=2)
  
}
ggarrange(plotlist=p.ls, common.legend=T, legend='right')

# Plot all CGIs together
cgi.all <- cgi.ls
for (i in 1:length(surveys)){
  cgi.all[[i]]$survey <- str_sub(surveys[i],1,4)
}
cgi.all <- do.call(rbind, cgi.all)
# Plot
p.base <- ggplot(data=dbs.Svmn.agg[[i]],
                 aes(x=x1/1000, y=x2/1000,fill=Kriging.Sv_mean.estim)) +
  geom_raster(show.legend=F) +
  scale_fill_gradientn(colours=c('white'),na.value="#A0A0A0") +
  labs(fill='Sv mean') +
  xlab("UTM Easting (km)") +
  ylab("UTM Northing (km)")
p1 <- p.base +
  geom_line(data=cgi.all[cgi.all$name == 'axis1',],
            aes(x=x/1000,y=y/1000, col=as.factor(survey)), inherit.aes=F, 
            linetype='solid', size=.35, alpha=.4) +
  geom_line(data=cgi.all[cgi.all$name == 'axis2',],
            aes(x=x/1000,y=y/1000, col=as.factor(survey)), inherit.aes=F, 
            linetype='solid', size=.35, alpha=.4) +
  geom_point(data=cgi.all[cgi.all$name == 'centre_gravity',],
            aes(x=x/1000,y=y/1000, col=as.factor(survey)), inherit.aes=F, 
            size=2.5) +
  labs(color='Year')
p2 <- p.base +
  geom_point(data=cgi.all[cgi.all$name == 'centre_gravity',],
             aes(x=x/1000,y=y/1000, col=as.factor(survey)), inherit.aes=F, 
             size=2.5) +
  labs(color='Year')
ggarrange(p2, p1, common.legend=T, legend='right')

###########################
# Plot penguin dive stats #
###########################
if (!file.exists('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')){
  dive.data <- tracks_diveStats(tracks.real, return_dive_df=T)
  saveRDS(dive.data, './data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
} else {
  dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
}

# Mean depth
p.d.mn <- ggplot(dive.data[[2]], aes(x=depth_mean, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Mean depth (m)') +
  ggtitle('Mean Dive Depth')

# Max depth
p.d.mx <- ggplot(dive.data[[2]], aes(x=depth_max, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Max depth (m)') +
  ggtitle('Max Dive Depth')

# Duration
p.d.dr <- ggplot(dive.data[[2]], aes(x=dive_duration, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Duration (s)') +
  ggtitle('Dive Duration')

# Vertical distribution again
# SV_mean
pSvmn <- ggplot(do.call(rbind, dbs.Svmn.vert), 
                aes(x=x3, y=Kriging.Sv_mean.estim, col=survey)) +
  geom_line() + geom_point() +
  scale_x_continuous(trans = "reverse") +
  coord_flip() +
  ggtitle('Sv mean')

# Plot
ggarrange(p.d.mn,  p.d.mx, p.d.dr, pSvmn, common.legend=T, legend='right')

# Repeat with better y scale
p.d.mn <- p.d.mn + ylim(0, .2)
p.d.mx <- p.d.mx + ylim(0, .12)
p.d.dr <- p.d.dr + ylim(0, .04)
# Plot
ggarrange(p.d.mn,  p.d.mx, p.d.dr, pSvmn, common.legend=T, legend='right')

#######################################################
# Now look at stats with dives longer than 20 seconds #
#######################################################
longCutoff <- 5
dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > longCutoff,]

# Mean depth
p.d.mn <- ggplot(dive.data.long, aes(x=depth_mean, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Mean depth (m)') +
  ggtitle('Mean Dive Depth')

# Max depth
p.d.mx <- ggplot(dive.data.long, aes(x=depth_max, col=factor(year(start)))) +
  geom_density() + 
  labs(col='Year') +
  xlab('Max depth (m)') +
  ggtitle('Max Dive Depth')

# Duration
p.d.dr <- ggplot(dive.data.long, aes(x=dive_duration, col=factor(year(start)))) +
  geom_density() +
  labs(col='Year') +
  xlab('Duration (s)') +
  ggtitle('Dive Duration')

ggarrange(p.d.mn,  p.d.mx, p.d.dr, pSvmn, common.legend=T, legend='right')

# Repeat with better y scale
p.d.mn <- p.d.mn + ylim(0, .15)
p.d.mx <- p.d.mx + ylim(0, .1)
p.d.dr <- p.d.dr + ylim(0, .04)
# Plot
ggarrange(p.d.mn,  p.d.mx, p.d.dr, pSvmn, common.legend=T, legend='right')

#####################################
# Plot dive density against Sv mean #
#####################################
# Plot side by side
Sv_vert <- do.call(rbind, dbs.Svmn.vert)
Sv_vert$year <- as.numeric(substr(Sv_vert$survey,1,4))

dive.data.long$year <- factor(year(dive.data.long$start)) 
p1 <- ggplot(dive.data.long, aes(x=depth_max, col=year)) +
  geom_density(size=1) + 
  labs(col='Year') +
  xlab('Depth (m)') + ylab('Density') +
  ggtitle('(a)') +
  coord_flip() +
  scale_x_reverse(limits=c(30,0))
  
p2 <- ggplot(Sv_vert, aes(x=x3, y=Kriging.Sv_mean.estim, color=as.factor(year))) +
  geom_line(size=1) +
  geom_point(size=2) +
  xlab(NULL) + ylab('SV mean') +
  ggtitle('(b)') +
  coord_flip() +
  scale_x_reverse(limits=c(30,0))
  
ggarrange(p1, p2, ncol=2, common.legend=T, legend='bottom')

# get Sv mean to deeper  (35m)
# Extract the top water of water column to depth
dbs.Svmn <- load.krig.dbs(var='Sv_mean', surveys=surveys, return_items=T)
dbs.Svmn <- lapply(dbs.Svmn, function(db) db[db$x3 <= 35 & db$Polygon,])
dbs.Svmn.vert <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=mean))
for (i in 1:length(surveys)){
  dbs.Svmn.vert[[i]]$survey <- surveys[i]
}
Sv_vert <- do.call(rbind, dbs.Svmn.vert)
Sv_vert$year <- as.numeric(substr(Sv_vert$survey,1,4))
row.names(Sv_vert) <- NULL

# Plot structure against dive desnity
# get desnity
p.build <- ggplot(dive.data.long, aes(x=depth_max, col=year)) +
  geom_density(size=1.3)
dive.density <- ggplot_build(p.build)
dive.density <- dive.density$data[[1]]
dive.density <- dive.density[,c('x','y','group')]
names(dive.density) <- c('depth','density','year')
dive.density$depth <- as.integer(round(dive.density$depth))
aggregate(dive.density, by=list(dive.density$year, dive.density$depth), mean)
j <- 2015
for (i in 1:5){
  dive.density$year[dive.density$year == i] <- j
  j <- j + 1
}
# aggregate density so we don't annoying "steps" in data
dive.density <- aggregate(dive.density, by=list(dive.density$year, dive.density$depth), 
                          mean)[,c(1,2,4)]
names(dive.density) <- c('year','depth','density')
# add in vertical depth
dive.density$Sv_mean <-
  unlist(mapply(function(d, s)
    Sv_vert$Kriging.Sv_mean.estim[d == Sv_vert$x3 & s == Sv_vert$year],
    dive.density$depth,
    dive.density$year))
dive.density$year <- factor(dive.density$year)
dive.density$Sv_linear <- Sv_mean.linear(dive.density$Sv_mean)

# plot
ggplot(dive.density, aes(x=density, y=Sv_linear, col=year)) +
  geom_point()
# smooth
ggplot(dive.density, aes(x=density, y=Sv_linear, col=year)) +
  geom_smooth()
# look at top 15m
ggplot(dive.density[dive.density$depth <= 15,], aes(x=density, y=Sv_linear, col=year)) +
  geom_smooth() + geom_point()

##############################################
# Look at horizontal vs vertical variability #
##############################################
# map the tracks and facet wrap for each year
tracks.map(dive.data[[1]], legend=F, facet_survey=T, colfactor = 'id', plot_line = F)
# now plot dive desnity by id for each year
ggplot(dive.data[[2]], aes(x=depth_max, col=id)) +
  geom_density(size=1.3) + 
  theme(legend.position = 'none') +
  xlab('Depth (m)') +
  ggtitle('Dive density') +
  coord_flip() +
  scale_x_reverse(limits=c(30,0)) +
  facet_wrap(~as.factor(year(start)))

#######################################################################
# GAM models for diving probability and vertical prey field structure #
#######################################################################

# Depth max.
dive.data.long$Sv_mean.dmx <- unlist(mapply(function(depth, year)
  Sv_vert$Kriging.Sv_mean.estim[Sv_vert$x3 == depth & Sv_vert$year == year],
  as.integer(round(dive.data.long$depth_max)),
  dive.data.long$year))
dive.data.long$Sv_linear.dmx <- Sv_mean.linear(dive.data.long$Sv_mean.dmx)

# Depth mean
dive.data.long$Sv_mean.dmn <- unlist(mapply(function(depth, year)
  Sv_vert$Kriging.Sv_mean.estim[Sv_vert$x3 == depth & Sv_vert$year == year],
  as.integer(round(dive.data.long$depth_mean)),
  dive.data.long$year))
dive.data.long$Sv_linear.dmn <- Sv_mean.linear(dive.data.long$Sv_mean.dmn)

# GAMS
gam.depthmean <- gam(depth_mean ~ s(Sv_linear.dmn, by=factor(year)) + year, 
                     data=dive.data.long)
summary(gam.depthmean)

plot(gam.depthmean, pages=1)


#################
# GAM attempt 2 #
#################
tracks.diveStats <- dive.data[[1]]
# add zeros
tracks.diveStats$diveDuration_sum[is.na(tracks.diveStats$diveDuration_sum)] <- 0
tracks.diveStats$diveDepth_mean[is.na(tracks.diveStats$diveDepth_mean)] <- 0
# Add SV_mean
Sv.proto <- mapply(function(depth, year) 
  Sv_vert$Kriging.Sv_mean.estim[Sv_vert$x3 == depth & Sv_vert$year == year],
  as.integer(round(tracks.diveStats$diveDepth_mean)),
  year(tracks.diveStats$dtUTC))
tracks.diveStats$Sv_mean.dmn <- unlist(lapply(Sv.proto, function (x) x[1]))
tracks.diveStats$Sv_linear.dmn <- Sv_mean.linear(tracks.diveStats$Sv_mean.dmn)
tracks.diveStats$year <- as.factor(year(tracks.diveStats$dtUTC))
# Drop missing positions
tracks.diveStats <- tracks.diveStats[complete.cases(tracks.diveStats[,c('lon','lat')]),]

# TEST GAM
for (s in str_sub(surveys,1,4)){
  mod <- gam(diveDuration_mean ~ s(lon, lat), data=tracks.diveStats[tracks.diveStats$year == s,]) 
  vis.gam(mod, plot.type='contour', main=s)
}



# for (depth in c(2,6,10,12,15,20,25)){
#   message(depth)
#   # Add Sv mean based on position
#   tracks.diveStats.set <- db.extract.covar.depth(tracks.diveStats, depth, convert2UTM=TRUE,
#                                              var='Sv_mean', method='slice')
#   tracks.diveStats.set$Sv_linear_krig <- Sv_mean.linear(tracks.diveStats.set$krig_Sv_mean)
#   # Build GAMS for dive effort for each depth...
#   gam <- gam(diveEffort ~ s(Sv_linear_krig, by=year) + year,
#                            data=tracks.diveStats.set)
#   plot(gam, pages=1, residuals=T, shade=T, pch=1, cex=1, main=depth)
# }

# Try with one dataframe
tracks.diveStats.set <- list()
i <- 1
for (depth in c(2,6,10,12,15,20,25)){
  message(depth)
  # Add Sv mean based on position
  tracks.diveStats.set[[i]] <- db.extract.covar.depth(tracks.diveStats, depth, convert2UTM=TRUE,
                                                 var='Sv_mean', method='slice')
  tracks.diveStats.set[[i]]$depth <- depth
  i <- i + 1
}
tracks.diveStats.set <- do.call(rbind, tracks.diveStats.set)
tracks.diveStats.set <- rows.and.levels(tracks.diveStats.set)
tracks.diveStats.set$diveEffort <- tracks.diveStats.set$diveDuration_mean*tracks.diveStats.set$diveCount
tracks.diveStats.set$Sv_linear_krig <- Sv_mean.linear(tracks.diveStats.set$krig_Sv_mean)
# Build GAM for dive effort for each depth...
for (s in str_sub(surveys,1,4)){
  mod <- gam(diveDuration_sum ~ te(Sv_linear_krig, depth),
             data=tracks.diveStats.set[tracks.diveStats.set$year == s,])
  vis.gam(mod, plot.type='contour', main=s, color='heat')
}



