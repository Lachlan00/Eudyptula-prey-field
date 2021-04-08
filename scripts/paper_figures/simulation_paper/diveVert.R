# Post hoc analysis to see why 2017 is different
library(ggplot2)
library(cmocean)
library(scales)
source('scripts/eudyptula.R')

# Paramaters
depth = 30 # Depth to extract Sv_mean data

# load databases
surveys <- c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
dbs.Svmn <- load.krig.dbs(var='Sv_mean', surveys=surveys, return_items=T)
# Extract the top water of water column to depth
dbs.Svmn <- lapply(dbs.Svmn, function(db) db[db$x3 <= depth & db$Polygon,])

# Side caclulation (report total Sv mean)
for (i in 1:length(dbs.Svmn)){
  message(names(dbs.Svmn)[i])
  print(mean(dbs.Svmn[[i]]$Kriging.Sv_mean.estim))
  print(sd(dbs.Svmn[[i]]$Kriging.Sv_mean.estim))
}

# Aggregate in 2 dimesnions
dbs.Svmn.vert <- lapply(dbs.Svmn, function(db) aggregate(db, by=list(db$x3), FUN=mean))

# Side calcs 
# Peak mean and depth for 2015
dbs.Svmn.vert$`2015_S1`[which.max(dbs.Svmn.vert$`2015_S1`$Kriging.Sv_mean.estim),]
# 2016 8-24 m depth range
mean(dbs.Svmn.vert$`2016_S2`[dbs.Svmn.vert$`2016_S2`$x3 > 7 & 
                             dbs.Svmn.vert$`2016_S2`$x3 <= 23,
                             'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2016_S2`[dbs.Svmn.vert$`2016_S2`$x3 > 7 & 
                               dbs.Svmn.vert$`2016_S2`$x3 <= 23,
                             'Kriging.Sv_mean.estim'])
# 2017 summaries
mean(dbs.Svmn.vert$`2017_S2`[dbs.Svmn.vert$`2017_S2`$x3 > 20, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2017_S2`[dbs.Svmn.vert$`2017_S2`$x3 > 20, 'Kriging.Sv_mean.estim'])
# 2018 summaries
mean(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 <= 5, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 <= 5, 'Kriging.Sv_mean.estim'])
mean(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 >= 20, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2018_S2`[dbs.Svmn.vert$`2018_S2`$x3 >= 20, 'Kriging.Sv_mean.estim'])
dbs.Svmn.vert$`2018_S2`[which.min(dbs.Svmn.vert$`2018_S2`$Kriging.Sv_mean.estim),]
dbs.Svmn.vert$`2018_S2`[which.max(dbs.Svmn.vert$`2018_S2`$Kriging.Sv_mean.estim),]
# 2019
mean(dbs.Svmn.vert$`2019_S2`[dbs.Svmn.vert$`2019_S2`$x3 <= 4, 'Kriging.Sv_mean.estim'])
sd(dbs.Svmn.vert$`2019_S2`[dbs.Svmn.vert$`2019_S2`$x3 <= 4, 'Kriging.Sv_mean.estim'])

# Add survey name to each database set
for (i in 1:length(surveys)){
  dbs.Svmn.vert[[i]]$survey <- surveys[i]
}

# load penguin dive data
dive.data <- readRDS('./data/analysis_datasets/cumsum/misc/posthoc_divedata.rds')
# dives longer than 5 seconds
longCutoff <- 5
dive.data.long <- dive.data[[2]][dive.data[[2]]$dive_duration > longCutoff,]
dive.data.long$year <- year(dive.data.long$start)

# Load real tracks to filter with
tracks.real <- readRDS('./data/analysis_datasets/cumsum/tracks_real.rds')
dive.data.long <- dive.data.long[as.character(dive.data.long$id) %in% 
                                 as.character(tracks.real$id),]

# Percentage dives < 30 m adn < 10 m
nrow(dive.data.long[dive.data.long$depth_max < 30,])/ nrow(dive.data.long)
nrow(dive.data.long[dive.data.long$depth_max < 10,])/ nrow(dive.data.long)

#####################################
# Plot dive density against Sv mean #
#####################################
# grey col
colGrey <- '#757575'
# Plot side by side
Sv_vert <- do.call(rbind, dbs.Svmn.vert)
Sv_vert$year <- as.numeric(substr(Sv_vert$survey,1,4))

# #--- VERT ---#
# dive.data.long$year <- factor(year(dive.data.long$start)) 
# p1 <- ggplot(dive.data.long, aes(x=depth_max)) +
#   geom_density(size=1, col=colGrey) + 
#   labs(col='Year') +
#   xlab('Depth (m)') + ylab('Density') +
#   ggtitle('a') +
#   coord_flip() +
#   scale_x_reverse(limits=c(30,0)) +
#   facet_wrap(~year, nrow=5) +
#   theme(legend.position = 'none') +
#   theme_bw() + grids(linetype = "dashed") +
#   theme(text = element_text(size=14))
# 
# p2 <- ggplot(Sv_vert, aes(x=x3, y=Kriging.Sv_mean.estim)) +
#   geom_line(size=1, col=colGrey) +
#   geom_point(size=2, col=colGrey) +
#   xlab(NULL) + ylab('Sv mean') +
#   ggtitle('b') +
#   coord_flip() +
#   scale_x_reverse(limits=c(30,0)) +
#   facet_wrap(~year, nrow=5) +
#   theme(legend.position = 'none') +
#   theme_bw() + grids(linetype = "dashed") +
#   theme(text = element_text(size=14))
# 
# p.vert <- ggarrange(p1, p2, ncol=2)
# ggsave(p.vert, filename='output/paper_figures/simulation_paper/diveVert_wrapVert.png',
#        dpi=300, width=1800/300, height=4000/300)

#--- SAME AXIS ---#
# try combine
p.vert <- ggplot(data=NULL) +
  geom_density(data=dive.data.long, mapping=aes(x=depth_max), 
               size=1, col=colGrey, alpha=.34, fill='grey') +
  geom_line(data=Sv_vert, mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17)))), 
            size=1, col='#c0392b', alpha=.8) +
  geom_point(data=Sv_vert, mapping=aes(x=x3, y=(rescale(Kriging.Sv_mean.estim, to=c(0, 0.17)))), 
            size=2, col='#c0392b', alpha=.8) +
  ylab('Density of dive depths') + xlab('Depth (m)') +
  scale_y_continuous(sec.axis = sec_axis(~rescale(.,
      to=c(min(Sv_vert$Kriging.Sv_mean.estim), max(Sv_vert$Kriging.Sv_mean.estim))), 
      name="Acoustic density (Sv)")) +
  coord_flip() +
  scale_x_reverse(limits=c(30,0)) +
  facet_wrap(~year, strip.position='right') +
  theme_bw() + grids(linetype = "dashed") +
  theme(text = element_text(size=14)) +
  theme(strip.background = element_rect(fill="#f2f2f2", colour="black"),
        strip.placement = "outside", 
        strip.text = element_text(size=12))

ggsave(p.vert, filename='output/paper_figures/simulation_paper/diveVert_sameAxis.png',
      dpi=300, width=3500/300, height=3000/300)

# RESCALE EXAMPLE used for secondary axis
# Sv_vert$Kriging.Sv_mean.estim
# rescale(rescale(Sv_vert$Kriging.Sv_mean.estim, to=c(0,0.16)), 
#         to=c(min(Sv_vert$Kriging.Sv_mean.estim), max(Sv_vert$Kriging.Sv_mean.estim)))


# #--- ORIG ---#
# p1 <- ggplot(dive.data.long, aes(x=depth_max, col=year)) +
#   geom_density(size=1) + 
#   labs(col='Year') +
#   xlab('Depth (m)') + ylab('Density') +
#   ggtitle('(a)') +
#   coord_flip() +
#   scale_x_reverse(limits=c(30,0)) +
#   theme_bw() + grids(linetype = "dashed") +
#   theme(text = element_text(size=14))
# 
# p2 <- ggplot(Sv_vert, aes(x=x3, y=Kriging.Sv_mean.estim, color=as.factor(year))) +
#   geom_line(size=1) +
#   geom_point(size=2) +
#   xlab(NULL) + ylab('SV mean') +
#   ggtitle('(b)') +
#   coord_flip() +
#   scale_x_reverse(limits=c(30,0)) +
#   theme_bw() + grids(linetype = "dashed") +
#   theme(text = element_text(size=14))
# 
# p.vert <- ggarrange(p1, p2, ncol=2, common.legend=T, legend='bottom')
# ggsave(p.vert, filename='output/paper_figures/simulation_paper/diveVert_orig.png',
#        dpi=300, width=3000/300, height=2500/300)
