# Check state sequences
library(momentuHMM)
library(reshape2)
source('scripts/Eudyptula.R')

# load HMM models
mod.HMM <- readRDS('./data/simulations/HMM_models/HMM.rds')
tracks <- mod.HMM$data
tracks$state <- viterbi(mod.HMM)
names(tracks)[c(1,4,5,6)] <- c('id','lon','lat','survey')
# Filter tracks that leave the survey area like in previous analysis
tracks <- inside.survey.zone(tracks, threshold=10, 
                             plot.map=F, plot.title='Real Tracks')
tracks$survey <- substr(tracks$survey,1,4)
tracks$state[tracks$state == 1] <- "Foraging"
tracks$state[tracks$state == 2] <- "Travel"
tracks$state <- as.factor(tracks$state)

# Calc proportion
track.ls <- split(tracks, tracks$id)
for (i in 1:length(track.ls)){
  dat <- track.ls[[i]]
  countState <- table(dat$state)/nrow(dat)
  row <- data.frame(id=dat$id[1], 
                    survey=dat$survey[1],
                    foraging=countState[1],
                    travel=countState[2])
  track.ls[[i]] <- row
}
track.data <- do.call(rbind, track.ls)
tracks.data.mlt <- melt(track.data)
tracks.data.mlt <- tracks.data.mlt[complete.cases(tracks.data.mlt),]

# Bar plot
ggplot(tracks.data.mlt, aes(x=survey, fill=variable, y=value)) +
  geom_boxplot(position="dodge") +
  labs(fill='State', x='Year', y='Proportion') +
  theme_bw() + grids(linetype = "dashed")

#####
# RLE
#####
# Calc rle
tracks.stateNA <- tracks
tracks.stateNA$state[is.na(tracks.stateNA$lon)] <- NA
tracks.stateNA <- rows.and.levels(tracks.stateNA)
track.ls <- split(tracks.stateNA, tracks.stateNA$id)
rle.ls <- list()
out2.ls <- list()
for (i in 1:length(track.ls)){
  year <- as.numeric(substr(track.ls[[i]]$survey,1,4))
  res <- rle(as.numeric(track.ls[[i]]$state))
  res <- data.frame(length=res$lengths, values=res$values)
  res <- res[complete.cases(res),]
  res <- res[res$values == 1,]
  res <- list(foraging.length=mean(res$length), rle=res, year=year[1])
  rle.ls[[i]] <- res
  out2.ls[[i]] <- data.frame(year=year[1], lengths=res$rle$length)
}
out2 <- do.call(rbind, out2.ls)
# mean rle
out <- data.frame(year=unlist(lapply(rle.ls, function(x) x$year)),
                  meanLen=unlist(lapply(rle.ls, function(x) x$foraging.length)))

# plot mean lengths for tracks
ggplot(out, aes(x=as.factor(year), y=meanLen)) +
  geom_boxplot()

# plot mean lengths for all bouts
ggplot(out2, aes(x=as.factor(year), y=lengths)) +
  geom_boxplot()


