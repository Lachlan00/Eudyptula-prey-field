library(readxl)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(gtools)
library(lazyWeave)
source('scripts/eudyptula.R')

# File list
dir <- './data/penguin/'
fn <- list.files(dir, pattern='*.xls*')

# Read in weight data from excel
peng.data <- lapply(fn, function(x) as.data.frame(read_excel(paste0(dir,x))))
peng.data <- lapply(peng.data, function(x)  x[3:nrow(x),])
peng.data <- data.frame(id=c(peng.data[[1]][,4],peng.data[[2]][,4],peng.data[[3]][,4]),
                        date=c(peng.data[[1]][,16],peng.data[[2]][,18],peng.data[[3]][,18]),
                        weight1=c(peng.data[[1]][,13],peng.data[[2]][,15],peng.data[[3]][,15]),
                        weight2=c(peng.data[[1]][,15],peng.data[[2]][,17],peng.data[[3]][,17]))
peng.data <- peng.data[complete.cases(peng.data),]
peng.data$date <- as.Date(as.numeric(as.character(peng.data$date)), origin="1899-12-30")
peng.data$weight1 <- as.numeric(as.character(peng.data$weight1))
peng.data$weight2 <- as.numeric(as.character(peng.data$weight2))
peng.data$id <- as.character(peng.data$id)

# Read in gemma weight data
peng.gemma <- readRDS("./data/misc/annotated_accel_deployments.rds")
peng.gemma <- data.frame(id=as.character(peng.gemma$likely_gps_track), 
                         date=as.Date(peng.gemma$Date),
                         weight1=peng.gemma$`Start weight`,
                         weight2=peng.gemma$`End weight`)

# Read in alternative (non acclrometer focused)  weight data
fn <- list.files(paste0(dir,'/gemma/'))
peng.gemma <- lapply(fn, function(x) as.data.frame(read_excel(paste0(dir,'gemma/',x))))
peng.gemma <- data.frame(id=as.character(c(paste0(peng.gemma[[1]]$Nest,peng.gemma[[1]]$Sex),
                                           paste0(peng.gemma[[2]]$Nest,peng.gemma[[2]]$Sex))), 
                         date=as.Date(c(peng.gemma[[1]]$Date, peng.gemma[[2]]$Date)),
                         weight1=c(peng.gemma[[1]]$`Start weight`,peng.gemma[[2]]$`Start weight`),
                         weight2=c(peng.gemma[[1]]$`End weight`,peng.gemma[[2]]$`End weight`))

# merge
peng.data <- rbind(peng.gemma, peng.data)
peng.data <- peng.data[complete.cases(peng.data$date) & complete.cases(peng.data$weight1),]

# Basic stats  test (with pseudoreplication not accounted for)
aov <- lm(weight1 ~ as.factor(year(date)), data=peng.data)
summary(aov)
# Looking good, 2017 significantly varies but need to make sure we only sample individuals once

# Sample individuals once
peng.data.norep <- peng.data
peng.data.norep$weight2  <- NULL
# Separate to L and G
L.pengs <- peng.data.norep[grepl('m|f|u',peng.data.norep$id),]
L.pengs$id <- substr(L.pengs$id,1,4)
G.pengs <- peng.data.norep[!grepl('m|f|u',peng.data.norep$id),]
G.pengs$id <- as.character(G.pengs$id)
# Take mean of duplicates
# L
L.pengs <- aggregate(L.pengs,  by=list(L.pengs$id), mean)
L.pengs <-  L.pengs[,c(1,3,4)]
colnames(L.pengs)[1] <- 'id'
# G
#G.pengs$id <- as.character(G.pengs$id)
G.pengs <- aggregate(G.pengs,  by=list(paste0(G.pengs$id,'_',year(G.pengs$date))), mean)
G.pengs <-  G.pengs[,c(1,3,4)]
colnames(G.pengs)[1] <- 'id'
# Merge
peng.data.norep <- rbind(G.pengs, L.pengs)

# separate sex
male <- sapply(peng.data.norep$id, function(id) grepl('m|M', id))
peng.data.norep$sex <- 'f'
peng.data.norep$sex[male] <- 'm'

# Also load and plot acoustics
surveys = c('2015_S1','2016_S2','2017_S2','2018_S2','2019_S2')
krig.dfs <- load.krig.dbs('Sv_mean', surveys, return_items = T)
for (i in 1:length(krig.dfs)){
  krig.dfs[[i]]$year <- (2015:2019)[i]
}
krig.dfs <- do.call(rbind, krig.dfs)
krig.dfs <- krig.dfs[krig.dfs$x3 <= 30,]
krig.dfs$linear <- Sv_mean.linear(krig.dfs$Kriging.Sv_mean.estim)

# Basic stats  test (with pseudoreplication accounted for)
# aov <- lm(weight1 ~ as.factor(year(date)), data=peng.data.norep)
# Filter penguins
if (FALSE){
  peng.data.norep <- split(peng.data.norep, year(peng.data.norep$date))
  # 2015
  peng.data.norep$`2015` <- peng.data.norep$`2015`[
    peng.data.norep$`2015`$date >= as.Date('2015-09-30') &
    peng.data.norep$`2015`$date <= as.Date('2015-10-08'),]
  # 2016
  peng.data.norep$`2016` <- peng.data.norep$`2016`[
    peng.data.norep$`2016`$date >= as.Date('2016-10-22') &
      peng.data.norep$`2016`$date <= as.Date('2016-10-31'),]
  # 2017
  peng.data.norep$`2017` <- peng.data.norep$`2017`[
    peng.data.norep$`2017`$date >= as.Date('2017-10-04') &
      peng.data.norep$`2017`$date <= as.Date('2017-10-13'),]
  # 2018
  peng.data.norep$`2018` <- peng.data.norep$`2018`[
    peng.data.norep$`2018`$date >= as.Date('2018-09-30') &
      peng.data.norep$`2018`$date <= as.Date('2018-10-07'),]
  # 2019
  peng.data.norep$`2019` <- peng.data.norep$`2019`[
    peng.data.norep$`2019`$date >= as.Date('2019-09-27') &
      peng.data.norep$`2019`$date <= as.Date('2019-10-09'),]
  # Rebind
  peng.data.norep <- do.call(rbind, peng.data.norep)
  row.names(peng.data.norep) <- NULL
}

# STATS
peng.data.norep$year <- as.numeric(year(peng.data.norep$date))
grand.mean <- mean(peng.data.norep$weight1)
peng.data.norep$diffs <- peng.data.norep$weight1 - grand.mean
p.list.weight  <- list()
i <- 1
message('weight')
for (YEAR in 2015:2019){
  message(YEAR)
  p.list.weight[[i]] <- t.test(peng.data.norep$diffs[peng.data.norep$year == YEAR])
  #p.list.weight[[i]] <- glm(diffs ~ year*sex, data=peng.data.norep)
  print(p.list.weight[[i]]$p.value)
  i <- i + 1
}
peng.data.norep$year <- as.factor(peng.data.norep$year)
summary(glm(diffs ~ year*sex, data=peng.data.norep))

##########################################################################################
# Add in core foraginf rnagw
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

##########################################################################################
        
# Report numbers
# Basic stats  test (with pseudoreplication accounted for)
grand.mean <- mean(krig.dfs$linear)
krig.dfs$diffs <- krig.dfs$linear - grand.mean
p.list.Sv.inside <- list()
p.list.Sv.outside <- list()
i <- 1
message('Sv')
for (YEAR in 2015:2019){
  message(YEAR)
  p.list.Sv.inside[[i]] <- wilcox.test(krig.dfs$diffs[krig.dfs$year == YEAR & krig.dfs$over],0)
  p.list.Sv.outside[[i]] <- wilcox.test(krig.dfs$diffs[krig.dfs$year == YEAR & !krig.dfs$over],0)
  message('inside ',p.list.Sv.inside[[i]]$p.value)
  message('outside ',p.list.Sv.outside[[i]]$p.value)
  i <- i + 1
}

summary(glm(diffs ~ year, data=krig.dfs))


# Plot based on year
text.weight <- pvalString(unlist(lapply(p.list.weight, function(x) x$p.value)))
text.weight <- data.frame(x=as.factor(2015:2019), y=-230, label=text.weight)
N.peng <- table(year(peng.data.norep$date))
N.peng <- data.frame(x=as.factor(2015:2019), 
                     y=-230, #c(0,30,-20,42,-1), 
                     label=as.vector(N.peng),
                     m=as.vector(table(year(peng.data.norep$date[peng.data.norep$sex == 'm']))),
                     f=as.vector(table(year(peng.data.norep$date[peng.data.norep$sex == 'f']))))
N.peng$label.sex <- paste0(N.peng$label,'\n(F',N.peng$f,':M',N.peng$m,')')

###########################################################
# Use updated code for new p1 plot (peng_wts.R) #
###########################################################
# p1 <- ggplot(peng.data.norep, aes(x=as.factor(year(date)), y=diffs, fill=as.factor(toupper(sex)))) +
#   geom_hline(yintercept = 0, color='darkred', linetype='dashed', size=.7) +
#   geom_boxplot( alpha=1) + #fill='#E59E7C'
#   labs(x=NULL, y='Weight anomaly (g)') +
#   theme_bw() + grids(linetype = "dashed") +
#   theme(text = element_text(size=14)) +
#   #theme(legend.position = "none") +
#   xlab(NULL) +
#   ylim(-240,230) +
#   labs(fill='Sex') +
#   # geom_label(data=text.weight, mapping=aes(x=x, y=y, label=label), 
#   #            inherit.aes = F, fill='#E59E7C', alpha=1) +
#   geom_text(data=N.peng, mapping=aes(x=x, y=y, label=paste0('N=',label.sex)), 
#              inherit.aes = F)

text.Sv.inside <- pvalString(unlist(lapply(p.list.Sv.inside, function(x) round(x$p.value,2))))
text.Sv.inside <- data.frame(x=as.factor(2015:2019), y=-3e-06, label=text.Sv.inside, over='Inside')
text.Sv.outside <- pvalString(unlist(lapply(p.list.Sv.outside, function(x) round(x$p.value,2))))
text.Sv.outside <- data.frame(x=as.factor(2015:2019), y=-3e-06, label=text.Sv.outside,  over='Outside')
text.Sv <- rbind(text.Sv.inside, text.Sv.outside)
text.Sv$x <- c(1:5 - .18, 1:5 + .18)


bool <- krig.dfs$over
krig.dfs$over[bool] <- 'Inside'
krig.dfs$over[!bool] <- 'Outside'

p2 <- ggplot(krig.dfs, aes(x=as.factor(year), y=diffs, fill=over)) +
  geom_hline(yintercept = 0, color='darkred', linetype='dashed', size=.7) +
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  labs(x=NULL, y='Sv mean anomaly (linear)', fill="Sex") +
  theme_bw() + grids(linetype = "dashed") +
  theme(text = element_text(size=18)) +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(-3.5e-6, 3.3e-5)) +
  #lims(y=c(-3e-6, 3e-6)) +
  geom_label(data=text.Sv, mapping=aes(x=x, y=y, label=label, fill=over), 
             inherit.aes = F, alpha=.6) +
  guides(label='none')

p3 <- ggarrange(p1, p2, ncol=1, common.legend = T)

ggsave(p3, filename='output/paper_figures/simulation_paper/penguinWeight.png',
       dpi=300, width=3750/300, height=3000/300, bg='white')



