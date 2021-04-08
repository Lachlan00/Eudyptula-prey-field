# This script merges all of Gemma and Lachlan's penguin track data
# and creates a pre processed clean data RDS file
library(beepr)
setwd("~/Development/PhD/repos/Eudyptula")
source("scripts/eudyptula.R")

#######################
# read and merge data #
#######################
dfL <- readRDS('data/gps/penguin_tracks_lachlan_RAW.rds')
dfG <- readRDS('data/gps/penguin_tracks_gemma_RAW.rds')
df <- do.call('rbind', list(dfG, dfL))
saveRDS(df, 'data/gps/penguin_tracks_all_RAW.rds')

##############
# clean data #
##############
# remove all points > 150 km from Montague
df <- range.clean(df, lon_orig=150.2269, lat_orig=-36.25201, km=150)
# trim tracks to foraging trip
df <- ocean_points_trim(df, return_ocean=TRUE)
# separate multiple foraging trips in tracks
df <- split_trips(df)
# remove rafting (points close to island)
# note this does not remove all rafting
df <- shapefile.filter(df)
# remove any remaining land points
df <- remove.land(df, calc_ocean=TRUE)
# clean duplicate point errors
df <- duplicate.error.clean(df, gap_threshold=300)
# Recursively filter out points with high velocites
df <- speed.clean(df, maxSpeed=20)

#############
# Save data #
#############
saveRDS(df, 'data/gps/penguin_tracks_all_CLEAN.rds')

###############################
# Make a regularised versions #
###############################
df_1m <- crawl_predict_all(df, 1)
saveRDS(df_1m, 'data/gps/penguin_tracks_all_REG-1m.rds')
df_5m <- crawl_predict_all(df, 5)
saveRDS(df_5m, 'data/gps/penguin_tracks_all_REG-5m.rds')
df_10m <- crawl_predict_all(df, 10)
saveRDS(df_10m, 'data/gps/penguin_tracks_all_REG-10m.rds')

beep()
