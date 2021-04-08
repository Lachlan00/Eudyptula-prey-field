# This is an alternative penguin cleaning pipeline that trims the tracks to land points
# so summary stats are better calcuated

# This script merges all of Gemma and Lachlan's penguin track data
# and creates a pre processed clean data RDS file
library(beepr)
source("scripts/eudyptula.R")

#######################
# read and merge data #
#######################
dfL <- readRDS('data/gps/penguin_tracks_lachlan_RAW.rds')
dfG <- readRDS('data/gps/penguin_tracks_gemma_RAW.rds')
df <- do.call('rbind', list(dfG, dfL))

##############
# clean data #
##############
# remove all points > 150 km from Montague
df <- range.clean(df, lon_orig=150.2269, lat_orig=-36.25201, km=150)
# trim tracks to foraging trip
df <- ocean_points_trim(df, return_ocean=TRUE, trim2land=TRUE)
# separate multiple foraging trips in tracks
df <- split_trips(df)
# clean duplicate point errors
df <- duplicate.error.clean(df, gap_threshold=300)
# Recursively filter out points with high velocites
df <- speed.clean(df, maxSpeed=20)

#############
# Save data #
#############
saveRDS(df, 'data/gps/penguin_tracks_all_CLEAN_sumStats.rds')

# ##############################
# # Make a regularised version #
# ##############################s
df_5m <- crawl_predict_all(df, 5)
saveRDS(df_5m, 'data/gps/penguin_tracks_all_REG-5m_sumStats.rds')

beep()
