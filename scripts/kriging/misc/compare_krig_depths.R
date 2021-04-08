# Compare krig field depths between surveys
library(ggplot2)
library(cmocean)
library(pbapply)
library(cmocean)

# depth max
max_depth <- 40

# Set surveys to test
surveys <- c('2015_S1',
             '2016_S1','2016_S2',
             '2017_S1','2017_S2',
             '2018_S1','2018_S2',
             '2019_S1','2019_S2')

# Load databases
dir <- './kriging/output/'
dbTemp <- pblapply(surveys, function(s) readRDS(paste0(dir,'CTD/models/3Dkrig_CTDdb_temperature_',
                                                       s,'.rds')))
dbSalt <- pblapply(surveys, function(s) readRDS(paste0(dir,'CTD/models/3Dkrig_CTDdb_salinity_',
                                                       s,'.rds')))
dbAggs <- pblapply(surveys, function(s) readRDS(paste0(dir,'agg/models/3Dkrig_aggdb_Sv_mean_',
                                                       s,'.rds')))

# remove non-ploygon
dbTemp <- pblapply(dbTemp, function(db) db@items[db@items$Polygon,])
dbSalt <- pblapply(dbSalt, function(db) db@items[db@items$Polygon,])
dbAggs <- pblapply(dbAggs, function(db) db@items[db@items$Polygon,])


# Aggregate models across depth
dbTemp <- pblapply(dbTemp, function(db) aggregate(db, 
                                                  by=list(db$x3),
                                                  FUN=function(x) mean(x, na.rm=T)))
dbSalt <- pblapply(dbSalt, function(db) aggregate(db, 
                                                  by=list(db$x3),
                                                  FUN=function(x) mean(x, na.rm=T)))
dbAggs <- pblapply(dbAggs, function(db) aggregate(db, 
                                                  by=list(db$x3),
                                                  FUN=function(x) mean(x, na.rm=T)))

# Attach survey info to models
for (i in 1:length(surveys)){
  dbTemp[[i]]$survey <- surveys[i]
  dbSalt[[i]]$survey <- surveys[i]
  dbAggs[[i]]$survey <- surveys[i]
}

# Merge data frames
dbTemp <- do.call(rbind, dbTemp)
dbSalt <- do.call(rbind, dbSalt)
dbAggs <- do.call(rbind, dbAggs)

# Set max depth
dbTemp <- dbTemp[dbTemp$x3 <= max_depth,]
dbSalt <- dbSalt[dbSalt$x3 <= max_depth,]
dbAggs <- dbAggs[dbAggs$x3 <= max_depth,]

# Plot the models
# Temp
ggplot(dbTemp, aes(x=Kriging.temperature.estim, y=x3, color=Kriging.temperature.estim)) +
  geom_point() +
  scale_color_gradientn(colors=cmocean('thermal')(100)) +
  scale_y_reverse() +
  labs(x='Temperature (celcius)', y='Depth (m)', col='Temperature') +
  facet_wrap(~survey)

# Salt
ggplot(dbSalt, aes(x=Kriging.salinity.estim, y=x3, color=Kriging.salinity.estim)) +
  geom_point() +
  scale_color_gradientn(colors=cmocean('haline')(100)) +
  scale_y_reverse() +
  labs(x='Salinity (PSU)', y='Depth (m)', col='Salinity') +
  facet_wrap(~survey)

# Aggs
ggplot(dbAggs, aes(x=Kriging.Sv_mean.estim, y=x3, color=Kriging.Sv_mean.estim)) +
  geom_point() +
  scale_color_gradientn(colors=rev(cmocean('deep')(100))) +
  scale_y_reverse() +
  labs(x='Sv mean', y='Depth (m)', col='Sv mean') +
  facet_wrap(~survey)



