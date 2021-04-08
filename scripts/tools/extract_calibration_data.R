# provide calibration data for Martin

# 1. Obtain CTD cast closest to calibrartion datetimes.
# 2. Provivde average sound speed, temperature and salinity between 1 - 12 m depth.

# libraries
library(stringr)

# source scritps
source('./scripts/CTD_functions.R')

#########################################################
# 1. Get calibration timestamps from calibrartion files #
#########################################################
# scan directory for calibration files
cal_dir <- '/Volumes/LP_MstrData/master-data/survey/acoustics/acoustics_sorted'
dir_ls <- list.dirs(cal_dir)
dir_ls <- dir_ls[grepl('*raw/calibration$', dir_ls)]
# remove 2018_S1 and 2019_S1 calibrartion was done in second survey
dir_ls <- dir_ls[!grepl('2018_S1|2019_S1', dir_ls)]
# get survey names
survey_names <- unlist(lapply(dir_ls, function(x) strsplit(x, '/')[[1]][8]))
# make a data frame
survey_df <- data.frame(survey_id=survey_names, cal_dir=dir_ls, cal_start=NA, cal_end=NA)
# capture list
df_ls <- list()
# get stime of calibration files
for (i in 1:nrow(survey_df)){
  # get cal files
  file_ls <- list.files(as.character(survey_df$cal_dir[i]), pattern="*.raw$")
  # filter krill
  file_ls <- file_ls[!grepl('kril_', file_ls, fixed=TRUE)]
  # make into dataframe
  df <- data.frame(survey_id=survey_df$survey_id[i],
                   cal_file=file_ls,
                   cal_time=NA)
  # for each calibration file extract datetimes
  for (j in 1:nrow(df)){
    # extract time
    sdate <- str_extract(df$cal_file[j], "D[0-9]{8}")
    stime <- str_extract(df$cal_file[j], "T[0-9]{6}")
    df$cal_time[j] <- as.POSIXct(paste0(sdate, stime), "D%Y%m%dT%H%M%S", tz='UTC')
  }
  df_ls[[i]] <- df
}
# join all
survey_calibrations <- do.call('rbind', df_ls)
survey_calibrations$cal_time <- as.POSIXct(survey_calibrations$cal_time, origin='1970-01-01',tz='UTC')

##############################################
# 2. Find nearest CTD casts and extract data #
##############################################
# read in CTD data
CTD_data <- cast.reader(location_from_time=FALSE, drop_non_transect=FALSE,
                        drop_duplicates=FALSE, drop_no_location=FALSE)
CTD_df <- CTD_data[[1]]
CTD_meta <- CTD_data[[2]]

# for each calibrartion find casts within 1.5 hours
survey_calibrations$cast_id <- NA
for (i in 1:nrow(survey_calibrations)){
  time_diff <- abs(difftime(CTD_meta$cast_time_UTC, survey_calibrations$cal_time[i], units='secs'))
  casts <- CTD_meta$file_name[which(time_diff < 1.5*3600)]
  if (length(casts) >= 1){
    survey_calibrations$cast_id[i] <- as.character(casts[length(casts)])
  } else {
    # if not found get as close as we can
    time_diff <- abs(difftime(CTD_meta$cast_time_UTC, survey_calibrations$cal_time[i], units='secs'))
    survey_calibrations$cast_id[i] <- as.character(CTD_meta$file_name[which.min(time_diff)])
  }
}

####################################################################################
# 3. Extract average sound speed, density, salinity and temperature between 1-12 m #
####################################################################################
survey_calibrations[,c('time_diff_hours','sound_speed','density','temperature','salinity')] <- NA
# loop and get the variables
for (i in 1:nrow(survey_calibrations)){
  # get cast
  cast <- CTD_df[as.character(CTD_df$id) == survey_calibrations$cast_id[i],]
  # round depth
  cast$depth <- round(cast$depth)
  # get depth range
  cast <- cast[cast$depth >= 1 & cast$depth <= 12,]
  # get variables
  survey_calibrations$time_diff_hours[i] <- round(difftime(survey_calibrations$cal_time[i],
                                                     CTD_meta$cast_time_UTC[CTD_meta$file_name == 
                                                                              as.character(cast$id[1])],
                                                     units='hours'),2)
  survey_calibrations$sound_speed[i] <- mean(cast$sound_velocity)
  survey_calibrations$density[i] <- mean(cast$density)
  survey_calibrations$temperature[i] <- mean(cast$temperature)
  survey_calibrations$salinity[i] <- mean(cast$salinity)
}
row.names(survey_calibrations) <- NULL
# save output
write.csv(survey_calibrations, './output/survey_calibration/survey_calibration.csv', row.names = FALSE)






