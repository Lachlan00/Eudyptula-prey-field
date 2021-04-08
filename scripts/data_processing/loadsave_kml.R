# load all kml files from source drive and save as quick reading RDS file
# function is recursive
setwd("~/Development/PhD/repos/Eudyptula")
source("scripts/data_functions.R")
load.save.tracks.source('/Volumes/LP_MstrData/master-data/penguins/tracks/', 
                        output_fn='penguin_tracks_lachlan_RAW2.rds')

# test
df <- readRDS('data/gps/loaded_KMLtracks.rds')
