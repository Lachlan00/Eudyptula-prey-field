# A set of functions for simulating penguin movement
source('./scripts/HMM/HMM_functions.R')
source('./scripts/utilities.R')
source('./scripts/visual.R')
source('./scripts/eudyptula.R')
source('./scripts/simulations/simData_custom.R')
# Libraries

#############################################################
# Create resampled track simulations from a mometuHMM model #
#############################################################
# Returns a dataframe of simulated tracks
simulate.HMM <- function(mod, nbTracks=20, spawnPoint=c(150.2180, -36.26215), plotMap=FALSE,
                         landShp='./assets/shp/NSW-coastline/edited/polygon/NSW-coastline_WGS84'){
  # convert Spawnpoint to UTM (km)
  spawnPoint <- lonlat2UTM(spawnPoint[1], spawnPoint[2])/1000
  spawnPoint <- c(spawnPoint[1], spawnPoint[2])
  names(spawnPoint) <- NULL
  spawnPoint <- unlist(spawnPoint)
  
  # load coast shapefile
  landShp <- readOGR(dsn=path.expand(str_sub(landShp,1,-(nchar(basename(landShp)))-1)), 
          layer=basename(landShp))
  # convert to UTM
  landShp <- spTransform(landShp, CRS("+proj=utm +zone=56H +ellps=WGS84"))
  
  # set observations per animal
  obsPerAnimal <- c(min(summary(mod$data$ID)), max(summary(mod$data$ID)))
  
  # Run simulation
  sim.df <- simData_custom(model=mod, nbAnimals=nbTracks, 
                           initialPosition=spawnPoint,
                           obsPerAnimal=obsPerAnimal,
                           landShp=landShp,
                           states=TRUE)
  names(sim.df) <- c('id','step','angle','lon','lat','state')
  # add state names
  sim.df$state <- factor(sim.df$state, labels=mod$stateNames)
  # Convert back to lonlat
  sim.df[,c('lon','lat')] <- sim.df[,c('lon','lat')]*1000
  sim.df <- UTM2lonlat.df(sim.df)
  # make a map
  if (plotMap)
    print(tracks.map(sim.df, zoom=11, legend=FALSE))
  return(sim.df)
}

##########################
# Load source HMM models #
##########################
# load HMM models
sim.load.HMM.models <- function(years=2015:2019, dir='./data/simulations/HMM_models/'){
  # load model for each year
  model.ls <- list()
  for (i in 1:length(years)){
    model.ls[[i]] <- readRDS(paste0(dir,'HMM_',years[i],'.rds'))
  }
  # names models
  names(model.ls) <- years
  return(model.ls)
}

##########################
# Check model parameters #
##########################
# Checks differences between model distributions for each HMM model
check.HMM.dists <- function(models){
  # Get params
  stepDist <- lapply(models, function(mod) mod$mle$step)
  angleDist <- lapply(models, function(mod) mod$mle$angle)
  # Put into dataframe formats
  stepDist <- data.frame(means = unlist(lapply(stepDist, function(dist) dist[1,])),
                   sds = unlist(lapply(stepDist, function(dist) dist[2,])),
                   year = rep(years, each=2),
                   state = rep(names(stepDist[[1]][1,]), length(years)))
  row.names(stepDist) <- NULL
  angleDist <- data.frame(means = unlist(lapply(angleDist, function(dist) dist[1,])),
                          concentrations = unlist(lapply(angleDist, function(dist) dist[2,])),
                          year = rep(years, each=2),
                          state = rep(names(angleDist[[1]][1,]), length(years)))
  row.names(angleDist) <- NULL
  # plot the results
  # Steps
  pStep <- ggplot(stepDist, aes(x=year, color=state)) +
    geom_errorbar(aes(ymax=means + sds, ymin=means - sds),
                  position="dodge", linetype = "dashed") +
    geom_point(aes(y=means), position=position_dodge(width=0.9), size=3) +
    ggtitle('HMM Step Gamma Mean and SD') +
    labs(x='Year', y='Step (km)', color='State')
  # Angle mean
  pAngleMean <- ggplot(angleDist, aes(x=year, color=state)) +
    geom_point(aes(y=means), position=position_dodge(width=0.9), size=3) +
    ggtitle('HMM Angle von Mises Means') +
    labs(x='Year', y='Mean Angle (radians)', color='State') + 
    ylim(c(-1,1)) +
    theme(legend.position='none')
  # Angle concentrations
  pAngleConc <- ggplot(angleDist, aes(x=year, fill=state)) +
    geom_bar(aes(y=concentrations), position='dodge', stat="identity") +
    ggtitle('HMM Angle Von Mises Concentrations') +
    labs(x='Year', y='Concentration around mean', fill='State') +
    theme(legend.position='none')
  # Angles
  pAngle <- ggarrange(pAngleMean, pAngleConc, nrow=2)
  # All together
  p <- ggarrange(pStep, pAngle, common.legend=T, legend='right')
  return(p)
}

#####################################################
# Report number of tracks by year for sims and real #
#####################################################
sim.report.tracksCount <- function(tracks){
  type.ls <- split(tracks, tracks$type)
  type.ls <- lapply(type.ls, function(x) split(x, x$year))
  df.real <- data.frame(year=names(type.ls$real),
                        count=unlist(lapply(type.ls$real, 
                                            function(df) length(unique(df$id)))))
  df.sim <- data.frame(year=names(type.ls$sim),
                       count=unlist(lapply(type.ls$sim,
                                           function(df) length(unique(df$id)))))
  message('Real tracks:')
  print(df.real)
  message('Simulated tracks:')
  print(df.sim)
}
