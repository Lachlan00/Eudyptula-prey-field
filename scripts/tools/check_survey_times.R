# scan aocustic aggregation data and report time ranges for each directory
input_dir <- './data/acoustics/'

# get files
file_ls <- list.files(input_dir, pattern="gps.csv$", recursive=TRUE)
survey_id <- substr(file_ls, 1, 7) 

# report data
for (i in 1:length(file_ls)){
  # get min and max times
  df <- read.csv(file.path(input_dir,file_ls[i]))
  dt <- as.POSIXct(paste(df$GPS_date, df$GPS_time), tz='UTC')
  dt_min <- min(dt)
  dt_max <- max(dt)
  # report
  message('')
  message(survey_id[[i]])
  message(dt_min,'\t',dt_max)
}
