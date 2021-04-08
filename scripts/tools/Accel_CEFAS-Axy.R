# Compare accelerotry data between CEFAS and Axy-Trek
library(ggplot2)
library(patchwork)
source('scripts/eudyptula.R')

# data directories
dir_CEFAS <- "/volumes/LP_MstrData/master-data/penguins/tracks/2017/dual_CEFAS/raw/"
dir_Axy <- "/volumes/LP_MstrData/master-data/penguins/tracks/2017/csv/"

# load in track codes with double deployments
ids <- list.files(dir_CEFAS, pattern='^[0-9]{3}[a-z]{1}[0-9]{2}$')

p_ls <- list()
# create the comparisons for each deployment
for (i in 1:length(ids)){
  id <- ids[i]
  message('Processing ',id)
  # find file
  fn <- list.files(paste0(dir_CEFAS,id,'/'), pattern='.CSV$')
  # load CEFAS
  message('\tReading CEFAS data..')
  CEFAS <- CEFAS.reader(paste0(dir_CEFAS,id,'/',fn))
  # load Axy
  message('\tReading AXY data..')
  # find file
  fn <- list.files(paste0(dir_Axy), pattern='.csv$')
  fn <- fn[grep(id, fn)]
  # check if multi file and then read in
  if (length(fn) == 1){
    # read in pressure
    AXY <- read.csv(paste0(dir_Axy, fn))
  } else {
    AXY <- do.call('rbind', lapply(paste0(dir_Axy, fn), read.csv))
  }
  message('\tCleaning data..')
  # clean AXY
  AXY$dtUTC <- as.POSIXct(paste(AXY$Date, AXY$Time), '%d/%m/%Y %H:%M:%OS', tz='UTC')
  AXY$id <- substr(AXY$TagID[1],1,6)
  AXY <- AXY[,c("id","dtUTC","X","Y","Z","Pressure","Temp....C.")]
  colnames(AXY)[6:7] <- c('pressure','temp')
  
  # remove group 0 from CEFAS
  CEFAS <- CEFAS[CEFAS$dive_id != 0,]

  # Filter all logs outside of Axy times
  Axy_time_range <- c(AXY$dtUTC[1], AXY$dtUTC[nrow(AXY)])
  fast_log_counts <- table(CEFAS$dive_id)
  CEFAS <- CEFAS[CEFAS$dtUTC > Axy_time_range[1] & CEFAS$dtUTC < Axy_time_range[2],]
  
  # set seed
  set.seed(420)
  # Select 4 fast logs to sample
  fast_log_counts <- table(CEFAS$dive_id)
  # only consider dives with more than 200 points
  fast_log_counts <- fast_log_counts[fast_log_counts > 200]
  # randomly choose 4
  fast_log <- as.numeric(names(sample(fast_log_counts, 1)))
  # make the plot
  CEFAS_sub <- CEFAS[CEFAS$dive_id == fast_log,]
  CEFAS_sub <- CEFAS_sub[1:200,]
  time_range <- c(CEFAS_sub$dtUTC[1], CEFAS_sub$dtUTC[nrow(CEFAS_sub)])
  AXY_sub <- AXY[AXY$dtUTC > time_range[1] & AXY$dtUTC < time_range[2],]
  # stack so can facet wrap
  AX.dat <- data.frame(Axis=rep(c('X','Y','Z'), each=nrow(AXY_sub)))
  AX.dat$dtUTC <- rep(AXY_sub$dtUTC, 3)
  AX.dat$mag <- c(AXY_sub$X,AXY_sub$Y,AXY_sub$Z)
  CEFAS.dat <- data.frame(Axis=rep(c('X','Y','Z'), each=nrow(CEFAS_sub)))
  CEFAS.dat$dtUTC <- rep(CEFAS_sub$dtUTC, 3)
  CEFAS.dat$mag <- c(CEFAS_sub$X, CEFAS_sub$Y, CEFAS_sub$Z)
  message('\tMaking plots..')
  # plot Axy  
  axy_plt <- ggplot(AX.dat, aes(x=dtUTC, y=mag, color=Axis)) +
    geom_line() +
    ggtitle(paste(id,'AXY')) +
    #facet_grid(rows=vars(Axis), scales='free_y') +
    #theme(legend.position = 'None') +
    labs(x='',y='')
  cefas_plt <- ggplot(CEFAS.dat, aes(x=dtUTC, y=mag, color=Axis)) +
    geom_line() +
    ggtitle(paste(id,'CEFAS')) +
    #facet_grid(rows=vars(Axis), scales='free_y') +
    #theme(legend.position = 'None') +
    labs(x='',y='')
  p_ls[[i]] <- (axy_plt / cefas_plt)
}

# make final plot
n <- length(p_ls)
ncol <- floor(sqrt(n))
do.call("grid.arrange", c(p_ls, ncol=ncol))




