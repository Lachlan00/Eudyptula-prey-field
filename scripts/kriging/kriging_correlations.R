# kriging correlations
library(ggplot2)
library(ggcorrplot)
library(cmocean)

# depth range
#depthRange <- c(0,20)
# set surveys to test
surveys <- c('2015_S1',
             '2016_S1','2016_S2',
             '2017_S1','2017_S2',
             '2018_S1','2018_S2',
             '2019_S1','2019_S2')

# load databases
dir <- './kriging/output/'
dbTemp <- pblapply(surveys, function(s) readRDS(paste0(dir,'CTD/models/3Dkrig_CTDdb_temperature_',
                                                       s,'.rds')))
dbSalt <- pblapply(surveys, function(s) readRDS(paste0(dir,'CTD/models/3Dkrig_CTDdb_salinity_',
                                                       s,'.rds')))
dbAggs <- pblapply(surveys, function(s) readRDS(paste0(dir,'agg/models/3Dkrig_aggdb_Sv_mean_',
                                                       s,'.rds')))

# calculate correlations between databases
dbdf <- data.frame(Polygon=dbTemp[[1]]@items$Polygon)
dbdf$x3 <- dbTemp[[1]]@items$x3
for (i in 1:length(surveys)){
  dbdf[,paste0('temp_',surveys[i])] <- dbTemp[[i]]@items$Kriging.temperature.estim
  dbdf[,paste0('salt_',surveys[i])] <- dbSalt[[i]]@items$Kriging.salinity.estim
  dbdf[,paste0('Sv_',surveys[i])] <- dbAggs[[i]]@items$Kriging.Sv_mean.estim
}
poly <- dbdf$Polygon
dbdf$Polygon <- NULL
dbdf <- dbdf[poly,]
#dbdf <- dbdf[dbdf$x3 >= depthRange[1] & dbdf$x3 <= depthRange[2],]
dbdf$x3 <- NULL
dbdf <- dbdf[complete.cases(dbdf),]
# correlations
corr_mat <- cor(dbdf)
p <- ggcorrplot(corr_mat, colors=c('#0b3487','#ffffff','#87130b'), lab=T, type='upper')
x1 <- c(.5,seq(3.5, (length(surveys)-2)*3+.5, 3),24.5)
x2 <- c(c(.5,seq(3.5, (length(surveys)-2)*3+.5, 3))+3,26.5)
y1 <- c(.5,seq(2.5, (length(surveys)-2)*3+.5, 3), 23.5)
y2 <- c(c(-.5,seq(2.5, (length(surveys)-2)*3+.5, 3))+3,26.5)
for (i in 1:length(x1)){
  p <- p + geom_rect(aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
              data=data.frame(x1=x1[i], x2=x2[i], y1=y1[i], y2=y2[i]), inherit.aes=F,
              size=1.6, fill="#00000000", color="black")  +
    geom_line(aes(x=x, y=y), linetype=2, color='black', size=1.6,
      data=data.frame(x=c(.5,x2[i]), y=y2[i]), inherit.aes=F) +
    geom_line(aes(x=x, y=y), linetype=2, color='black', size=1.6,
              data=data.frame(x=x2[i], y=c(.5,y2[i])), inherit.aes=F)
}
print(p)

# Plot 2019_S1 by itself
corr_solo <- corr_mat[grepl('2019_S1', rownames(corr_mat)),
                      grepl('2019_S1', colnames(corr_mat))]
rownames(corr_solo) <- c('Temperature', 'Salinity', 'Sv mean')
colnames(corr_solo) <- c('Temperature', 'Salinity', 'Sv mean')
ggcorrplot(corr_solo, colors=cmocean('balance')(5)[2:4], lab=T, type='upper')
