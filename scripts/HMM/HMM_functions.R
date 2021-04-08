# functions for producing Hidden Markov Models
# library(moveHMM)
library(momentuHMM)
library(lubridate)
library(RANN)

# NOTE: moveHMM is not loaded directly to avoid function masking issues. 
# Call moveHMM functions with "::" (e.g. moveHMM::prepData(regdata)).

# load source scritps
source('./scripts/utilities.R')
source('./scripts/eudyptula.R')

################################
# momentuHMM  prepdata wrapper #
################################
# fit a HMM model using momentuHMM
prepData.wrapper <- function(regdata, coord_type='LL', coordNames=c('lon','lat')){
  # compute steps and angles
  colnames(regdata)[colnames(regdata) == 'id'] <- 'ID'
  data <- momentuHMM::prepData(regdata, type=coord_type, coordNames=coordNames)
  return(data)
}

############################
# Fit HMM with momentuHMM  #
############################
# fit a HMM model using momentuHMM
fit.momentuHMM <- function(regdata, nbStates, dist, Par0, formula=~1, standardiseCovs=TRUE,
                           stateNames=NULL, coord_type='LL', coordNames=c('lon','lat'),
                           fixPar=NULL, centers=NULL, filter.gaps=FALSE, estAngleMean=NULL,
                           hack0steps=TRUE, keepCols=NULL){
  message('Preparing data..')
  # get covariate names from formula
  covNames <- strsplit(as.character(formula)[2], ' + ', fixed=T)[[1]]
  if (covNames == '1')
    covNames <- NULL
  # compute steps and angles
  colnames(regdata)[colnames(regdata) == 'id'] <- 'ID'
  data <- momentuHMM::prepData(regdata, type=coord_type, coordNames=coordNames, 
                               covNames=covNames, centers=centers)
  # Hack step
  # any zero values for steps get  a very small value added
  if (hack0steps){
    if (sum(data$step == 0 & !is.na(data$step)) >=1){
      data$step[data$step == 0 & !is.na(data$step)] <- 
        data$step[data$step == 0 & !is.na(data$step)] + 1e-4
      warning('1e-4 added to zero steps to avoid zero mass error messsage')
    }
  }
  # filter gaps
  if (filter.gaps){
    data <- tracks.reg.filter.gaps(data, coordNames=c('x','y'))
  }
  # remove unwanted streams
  if (length(names(dist)) > 2){
    data_streams <- names(dist)[3:length(names(dist))]
  } else {
    data_streams <- NULL
  }
  colony_streams <- NULL
  if (!is.null(centers)){
    colony_streams <- paste0(row.names(centers),'.',c('dist','angle'))
  }
  data <- data[,c('ID','step','angle',data_streams,'x','y',colony_streams,covNames,keepCols)]
  # standardise covariate data more than one covariate
  if (as.character(formula)[2] != '1' & length(covNames) > 1 & standardiseCovs){
    for (var in covNames)
      data[,var] <- (data[,var] - mean(data[,var])) / sd(data[,var])
  }
  # fit 2-state model
  message('Fitting Hidden Markov Model..')
  m <- momentuHMM::fitHMM(data=data,
                          nbStates=nbStates,
                          dist=dist,
                          Par0=Par0, 
                          formula=formula,
                          stateNames=stateNames,
                          estAngleMean=estAngleMean,
                          fixPar=fixPar)
  # return model
  return(m)
}

#########################
# Fit HMM with moveHMM  #
#########################
# fit a HMM model using moveHMM
fit.moveHMM <- function(regdata, nbStates, stepPar0, anglePar0, formula=~1){
  message('Preparing data..')
  # get covariate names from formula
  covNames <- strsplit(as.character(formula)[2], ' + ', fixed=T)[[1]]
  # compute steps and angles
  colnames(regdata)[colnames(regdata) == 'id'] <- 'ID'
  colnames(regdata)[colnames(regdata) == 'lon'] <- 'x'
  colnames(regdata)[colnames(regdata) == 'lat'] <- 'y'
  data <- moveHMM::prepData(regdata)
  # standardise covariate data if needed
  if (as.character(formula)[2] != '1' & length(covNames) > 1){
    for (var in covNames)
      data[,var] <- (data[,var] - mean(data[,var])) / sd(data[,var])
  }
  # fit 2-state model
  message('Fitting Hidden Markov Model..')
  m <- moveHMM::fitHMM(data, nbStates=nbStates, stepPar0=stepPar0, 
              anglePar0=anglePar0, verbose=0, formula=formula)
  # decode states
  return(m)
}

###################################################
# Get covariate data from kriging model databases #
###################################################
# uses KD trees to find penguin location in kriging model databases and
# extract covariate data
tracks.attach.krig.covar <- function(tracks, covars=c('temperature','salinity','Sv_mean','indicator',
                                                      'bathymetry'),
                                     depth_colname='depth_percen85',
                                     db_dirs=c('./kriging/output/CTD/models/','./kriging/output/agg/models/',
                                               './kriging/output/indicator/models/'),
                                     bathy_db_fn=c('./kriging/output/bathymetry/models/'),
                                     survey_time_fn='./data/transects/survey_times.csv',
                                     extract_quantiles=FALSE, # deprecatred
                                     acoustic_multi_depth_extract=FALSE,
                                     depth_vector=c(5,10,15,20,30,40),
                                     calc_polygon=FALSE, # deprecated?
                                     verbose=FALSE,
                                     convertUTM=TRUE){
  # create UTM coordinates (filter for 2015 before function)
  if (convertUTM){
    tracks[,c('UTM_easting','UTM_northing')] <- lonlat2UTM.track(tracks)[,c('lon','lat')]
  } else {
    tracks[,c('UTM_easting','UTM_northing')] <- tracks[,c('lon','lat')]
  }
  # make rounded depth
  tracks$rdepth <- round(tracks[,depth_colname])
  # load survey times
  survey_times <- read.csv(survey_time_fn)
  survey_times$start_UTC <- as.POSIXct(survey_times$start_UTC, tz='UTC')
  survey_times$end_UTC <- as.POSIXct(survey_times$end_UTC, tz='UTC')
  # stack columns
  survey_times_2 <- survey_times
  survey_times_2$start_UTC <- survey_times_2$end_UTC
  survey_times <- rbind(survey_times, survey_times_2)
  survey_times$end_UTC <- NULL
  colnames(survey_times)[3] <- 'dtUTC'
  # find which survey each track should pull data from
  track_survey_links <- data.frame(track_id=unique(tracks$id))
  track_survey_links$mean_dtUTC <- as.POSIXct(sapply(track_survey_links$track_id, 
                                          function(id) mean(tracks$dtUTC[tracks$id == id])),
                                          origin='1970-01-01', tz='UTC')
  # find survey from time
  track_survey_links$survey_id <- sapply(track_survey_links$mean_dtUTC,
                                         function(mean_dt) survey_times$id[which.min(abs(
                                           difftime(survey_times$dtUTC, mean_dt, units='hour')))])
  # add linking file to tracks
  # quicker to use rle
  message('Linking tracks to survey files..')
  track_ids_rle <- rle(as.vector(tracks$id))
  survey_ids <- as.character(sapply(track_ids_rle$values, function(x) 
    track_survey_links$survey_id[track_survey_links$track_id == x]))
  tracks$survey_id <- rep(survey_ids, times=track_ids_rle$lengths)
  # get database file names
  file_ls <- do.call(c, lapply(db_dirs, list.files))
  # filter for data bases only 
  # (this step is obsolete now as I no longer store dataframes)
  file_ls <- file_ls[substr(file_ls, 11,12) == 'db' | substr(file_ls, 6,7) == 'db']
  # load the covariate data
  
  # NOTE: It is faster to avoid looping the variables as each model grid is the same 
  # Consequently indexes only need to be obtained once for each track.
  # UPDATE: kd-trees are so fast that we should just cycle variables anyway.
  survey_set <- unique(tracks$survey_id)
  track_ls <- split(tracks, tracks$survey_id)
  
  # check if bathymetry is in covars and if so remove and make separate section for attaching
  if ('bathymetry' %in% covars){
    add_bathy <- TRUE
    covars <- covars[!(covars=='bathymetry')]
    # load bathymetry database
    bathy_db <- readRDS(paste0(bathy_db_fn, 
                               list.files(bathy_db_fn, pattern='krig_bathymetry-db.rds')))@items
  } else {
    add_bathy=FALSE
  }
  
  # list for storing a database from each var (for mater making polygon filter)
  db_polysets <- list()
  # Get covariates for each variable/survey model pair
  for (var in covars){
    message('\nProcessing ',var)
    # load the kriging models
    # get subset of the file list that is the variable
    file_idx <- which(grepl(var, file_ls))
    # load all the models for the variable
    # TEMP HACK FOR ACOUSTICS
    dir_idx <- 1
    if (var == 'Sv_mean')
      dir_idx <- 2
    if (var == 'indicator')
      dir_idx <- 3
    message(paste('Loading',var,'models..'))
    db_ls <- pblapply(paste0(db_dirs[dir_idx],file_ls[file_idx]), function(x) readRDS(x)@items)
    # get a sample for polygons
    db_polysets[[which(var == covars)]] <- db_ls[[1]]
    # get survey names
    if (dir_idx == 1){
      survey_names <- lapply(file_ls[file_idx], function(s) strsplit(s,'_')[[1]][4:5])
    } else if (dir_idx == 2) {
      survey_names <- lapply(file_ls[file_idx], function(s) strsplit(s,'_')[[1]][5:6])
    } else if (dir_idx == 3){
      survey_names <- lapply(file_ls[file_idx], function(s) strsplit(s,'_')[[1]][3:4])
    }
    survey_names <- unlist(lapply(survey_names, function(s) substr(paste(s[1],s[2],sep='_'),1,7)))
      
    # go through survey set and get the correct database
    for (i in 1:length(survey_set)){
      survey <- survey_set[i]
      message('Processing ',survey)
      track_df <- track_ls[[i]]
      track_df$id <- droplevels(track_df$id)
      # cut missing values (in depth) from track
      # we only cut from the tail or head to avoid messing up the regularisation
      # there should be no missing values within the track if previous scripts have done there job.
      track_df <- split(track_df, track_df$id)
      for (j in 1:length(track_df)){
        depthna_rle <- rle(is.na(track_df[[j]]$rdepth))
        # cut tail
        if (depthna_rle$values[length(depthna_rle$values)]){
          track_df[[j]] <- track_df[[j]][1:(sum(depthna_rle$lengths) - 
                                                 depthna_rle$lengths[length(depthna_rle$lengths)]),]
          if (verbose){
            print('Tail trimmed!')
            print(as.character(track_df[[j]]$id[1]))
          }
        }
        # cut head (don't think this will ever need to happen)
        if (depthna_rle$values[1]){
          track_df[[j]] <- track_df[[j]][(depthna_rle$lengths[1]+1):nrow(track_df[[j]]),]
          if (verbose){
            print('Head trimmed!')
            print(as.character(track_df[[j]]$id[1]))
          }
        }
      }
      track_df <- do.call('rbind', track_df)
      track_df$id <- droplevels(track_df$id)
      row.names(track_df) <- NULL
      # load database
      db <- db_ls[[which(survey_names == survey)]]
      db <- db[db$Polygon,]
      row.names(db) <- NULL
      # make sure we do not query NA values (nn2 will fail)
      na.idx <- which(is.na(track_df$rdepth))
      track_df_nafiltered <- track_df[!is.na(track_df$rdepth),]
      # use a kd-tree to quickly find the index for the databases
      db_idx <- as.vector(nn2(db[,c('x1','x2','x3')], 
                              track_df_nafiltered[,c('UTM_easting', 'UTM_northing','rdepth')], k=1)$nn.idx)
      # sub in NAs back in to vector
      for (ii in na.idx)
        db_idx <- c(db_idx[1:(ii-1)],NA,db_idx[ii:length(db_idx)])
      # attach to track data_frame
      if (var != 'indicator'){
        track_df[,c(paste0('krig_',var))] <- db[db_idx, paste0('Kriging.',var,'.estim')]
      } else {
        classes <- c('zero','low','medium','large','extrem','maxi','class')
        track_df[,c(paste0('krig_indicator_',classes))] <- db[db_idx, classes]
      }
      # if (extract_quantiles & var == 'Sv_mean')
      #   track_df[,c(paste0('krig_',var,'_quant'))] <- db[db_idx, paste0('quants')]
      # if option is on also attach variables at multiple depth values
      if (acoustic_multi_depth_extract & var %in% c('Sv_mean','indicator')){
        for (depth in depth_vector){
          db_idx <- as.vector(nn2(db[db$x3 == depth,c('x1','x2')], 
                                  track_df[,c('UTM_easting', 'UTM_northing')], k=1)$nn.idx)
          # attach to track data_frame
          if (var != 'indicator'){
            track_df[,c(paste0('krig_',var,'_',depth))] <- db[db$x3 == depth,][db_idx, paste0('Kriging.',var,'.estim')]
          } else {
            track_df[,c(paste0('krig_indicator_',classes,'_',depth))] <- db[db$x3 == depth,][db_idx, classes]
          }
          # if (extract_qunatiles)
          #   track_df[,c(paste0('krig_',var,'_quant_',depth))] <- db[db$x3 == depth,][db_idx,'quants']
        }
      }
      # add to track list
      track_ls[[i]] <- track_df
    }
  }
  # merge tracks
  tracks <- do.call('rbind', track_ls)
  # and finally add bathymetry
  message('Adding bathymetry data to tracks..')
  if (add_bathy){
    bathy_idx <- as.vector(nn2(bathy_db[,c('x1','x2')], 
                     tracks[,c('UTM_easting', 'UTM_northing')], k=1)$nn.idx)
    tracks$bathymetry <- bathy_db$Kriging.depth.estim[bathy_idx]
    db_polysets[[length(db_polysets)+1]] <- bathy_db
  }
  if (calc_polygon){
    message('Calculating unified kriging field polygon..')
    # add polygon data (for later filtering)
    poly_df <- data.frame(poly1=rep(NA,nrow(tracks)))
    for (i in 1:length(db_polysets)){
      # check if x3 present
      if ('x3' %in% names(db_polysets[[i]])){
        poly_idx <- as.vector(nn2(db_polysets[[i]][,c('x1','x2','x3')],
                              tracks[,c('UTM_easting', 'UTM_northing','rdepth')], k=1)$nn.idx)
      } else {
        poly_idx <- as.vector(nn2(db_polysets[[i]][,c('x1','x2')],
                                  tracks[,c('UTM_easting', 'UTM_northing')], k=1)$nn.idx)
      }
      poly_df[,paste0('poly',i)] <- db_polysets[[i]][poly_idx, 'Polygon']
    }
    tracks$Polygon <- Reduce("&",poly_df)
  }
  # clean up
  tracks[,c('UTM_easting', 'UTM_northing','rdepth')] <- NULL
  colnames(tracks)[colnames(tracks) == 'survey_id'] <- 'krig_survey_id'
  row.names(tracks) <- NULL
  tracks$id <- droplevels(tracks$id)
  return(tracks)
}

##################################
# Remove rafting from track tail #
##################################
# Cuts out rafting period when penguins are about to return to colony
# Maintains time regularisation of tracks
# radius in km
# window in minutes
# Uses HMM model settings to isolate rafting periods and not transit periods
remove.rafting <- function(tracks, nbstates=2, stepPar0, anglePar0, 
                           center=c(150.2269, -36.25201), radius=2.5, 
                           window=60, plot_process=FALSE, verbose=FALSE){
  # get track time interval and calculate how many points needed to confirm window
  td <- as.numeric(difftime(tracks$dtUTC[2], tracks$dtUTC[1], units='mins'))
  window_adjusted <- floor(window/td)
  # make the cirle
  circle <- circlemaker(center, radius)
  # find which points are in radius
  tracks$in_radius <- as.logical(point.in.polygon(tracks$lon, tracks$lat, 
                                                  circle@coords[,1], circle@coords[,2]))
  # plot if specified
  if (plot_process){
    # load basemap
    basemap <- load.basemap('montague', 12) 
    # plot basemap and circle
    print(basemap +
            geom_point(data=tracks, aes(lon, lat, color=in_radius), alpha=0.3) +
            geom_polygon(data=fortify(circle), aes(long, lat), color='red', fill=NA) +
            ggtitle('Rafting Radius')+
            labs(color='In Radius'))
  }
  # work on each track independently
  tracks <- split(tracks, tracks$id)
  message('Removing rafting points...')
  # if not verbose use progress bar
  if (!verbose)
    pb <- txtProgressBar(min=0, max=length(tracks), style = 3)
  # loop through and trim rafting tails
  for (i in 1:length(tracks)){
    if (!verbose)
      setTxtProgressBar(pb, i)
    track <- tracks[[i]]
    row.names(track) <- NULL
    # run length encode in_radius
    track_rle <- rle(track$in_radius)
    # only consider track if it finishes inside radius
    if (!track_rle$values[length(track_rle$values)]){
      if (verbose)
        message('\t',track$id[1],' - 0 rafting points removed')
      next()
    }
    # if track finishes in radius consider the last portion inside the radius
    track_radius <- track[(nrow(track)-track_rle$lengths[length(track_rle$lengths)]+1):nrow(track),]
    # only consider track if it has more than 3 points in the radius!
    if (nrow(track_radius) <= 3){
      if (verbose)
        message('\t',track$id[1],' - 0 rafting points removed')
      next()
    }
    # perform an HMM model on radius section
    r.HMM <- tryCatch(suppressWarnings(
      suppressMessages(
        fit.moveHMM(track_radius[,c('id','lon','lat')], nbStates, stepPar0, anglePar0))),
      error=function(e) return(NA))
    # if HMM process failed abort
    if (!is(r.HMM,'moveHMM')){
      if (verbose)
        message('\t',track$id[1],' - 0 rafting points removed')
      next()
    }
    # attach states to track_radius
    states <- suppressWarnings(tryCatch(as.character(moveHMM::viterbi(r.HMM)), error=function(e) return(NA)))
    # if failed to decode states assume all are travel and move on
    # non robust way to check for failure... but it works in this context...
    if (length(states) == 1){
      if (verbose)
        message('\t',track$id[1],' - 0 rafting points removed')
      next()
    }
    # rename states
    states[states == "1"] <- "ARS"
    states[states == "2"] <- "Travel"
    states <- as.factor(states)
    # create a check vector to test if penguin is both in radius AND in ARS state. 
    check_vector <- states == 'ARS'
    check_vector <- c(rep(FALSE, nrow(track) - nrow(track_radius)), check_vector)
    # if all FALSE
    if (sum(check_vector) == 0){
      if (verbose)
        message('\t',track$id[1],' - 0 rafting points removed')
      next()
    }
    # rle the vector
    check_vector_rle <- data.frame(unclass(rle(check_vector)))
    # find first (starting from end) period that is greater than the window size and FALSE
    for (j in 1:nrow(check_vector_rle)){
      row <- check_vector_rle[seq(dim(check_vector_rle)[1],1),][j,]
      if (!row$values & row$lengths > window_adjusted){
        idx_cut <- sum(check_vector_rle$lengths[1:(nrow(check_vector_rle) - (j-1))]) + 1
        break()
      }
    }
    # filter the track
    if (plot_process){
      track$cut <- TRUE
      track[1:idx_cut,'cut'] <- FALSE
      print(suppressWarnings(suppressMessages(tracks.map(track, colfactor='cut'))))
      track$cut <- NULL
    }
    before_count <- nrow(track)
    tracks[[i]] <- track[1:idx_cut,]
    if (verbose)
      message('\t',track$id[1],' - ',before_count - nrow(tracks[[i]]),' rafting points removed')
  }
  if (!verbose)
    close(pb)
  # merge filtered tracks
  tracks <- do.call('rbind', tracks)
  row.names(tracks) <- NULL
  tracks$in_radius <- NULL
  return(tracks)
}
  
                                          



