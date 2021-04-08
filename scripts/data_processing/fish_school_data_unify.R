# Merge 2016_S2 and 2018_S1 aggegation, depth and GPS files into a single files 
# to unify processing across surveys
# Do each section separately as they are quite different in structure
# main directory
data_dir <- '/Volumes/LP_MstrData/master-data/survey/acoustics/acoustics_processed/'

###########
# 2016_S2 #
###########
# get file list
file.ls <- list.files(paste0(data_dir,'2016_S2'))
# separate into each type
agg.ls <- file.ls[grepl('^agg*', file.ls)]
GPS.ls <- file.ls[grepl('^GPS*', file.ls)]
depth.ls <- file.ls[grepl('^Seabed*', file.ls)]

# Aggregations
# load and merge all files (non ASIC characters in column names hence check.names)
schools <- do.call('rbind', lapply(paste0(data_dir, '2016_S2/', agg.ls), read.csv, check.names=F))
colnames(schools)[1] <- 'Region_ID'
schools_dt <- as.POSIXct(paste0(schools$Date_S, schools$Time_S), format='%Y%m%d %H:%M:%OS', tz='UTC')
schools <- schools[order(schools_dt),]
write.csv(schools, paste0(data_dir,'2016_S2/2016_S2_schools.csv'), row.names=F)

# GPS
# load and merge all files
GPS <- do.call('rbind', lapply(paste0(data_dir, '2016_S2/', GPS.ls), read.csv))
GPS_dt <- as.POSIXct(paste0(GPS$GPS_date, GPS$GPS_time,'.',GPS$GPS_milliseconds), format='%Y-%m-%d %H:%M:%OS', tz='UTC')
GPS <- GPS[order(GPS_dt),]
write.csv(GPS, paste0(data_dir,'2016_S2/2016_S2.gps.csv'), row.names=F)

# Depth
# load and merge all files
depth <- do.call('rbind', lapply(paste0(data_dir, '2016_S2/', depth.ls), read.csv))
depth_dt <- as.POSIXct(paste0(depth$Ping_date, depth$Ping_time,'.',depth$Ping_milliseconds), format='%Y-%m-%d %H:%M:%OS', tz='UTC')
depth <- depth[order(depth_dt),]
write.csv(depth, paste0(data_dir,'2016_S2/2016_S2_Editable seabed.depth.csv'), row.names=F)

###########
# 2018_S1 #
###########
# get file list
file.ls <- list.files(paste0(data_dir,'2018_S1'))
# depth.ls <- file.ls[grepl('^Seabed*', file.ls)]

# Aggregations
# In this case the file was previously merged so we only need to load "20180822_AGG.csv"
schools <- read.csv(paste0(data_dir, '2018_S1/20180822_AGG.csv'))
schools_dt <- as.POSIXct(paste0(schools$Date_S, schools$Time_S), format='%Y%m%d %H:%M:%OS', tz='UTC')
schools <- schools[order(schools_dt),]
write.csv(schools, paste0(data_dir,'2018_S1/2018_S1_schools.csv'), row.names=F)

# GPS
# In this case the GPS is attached to the school aggregations so we can skip this

# Depth
# In this case the file was previously merged so we only need to load "20180822_seabed.csv"
depth <- read.csv(paste0(data_dir, '2018_S1/20180822_seabed.csv'))
depth_dt <- as.POSIXct(paste0(depth$Ping_date, depth$Ping_time,'.',depth$Ping_milliseconds), format='%Y-%m-%d %H:%M:%OS', tz='UTC')
depth <- depth[order(depth_dt),]
write.csv(depth, paste0(data_dir,'2018_S1/2018_S1_Editable seabed.depth.csv'), row.names=F)

