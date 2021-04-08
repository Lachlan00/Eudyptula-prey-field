# Process Gemma's tracks
setwd("~/Development/PhD/repos/Eudyptula")
source("scripts/eudyptula.R")
input_dir = '/Volumes/LP_MstrData/master-data/gemma/all_penguin_tracks/raw/'
output_dir = '/Volumes/LP_MstrData/master-data/gemma/all_penguin_tracks/processed/penguin_tracks_gemma_RAW.rds'
meta_out = '/Volumes/LP_MstrData/master-data/gemma/all_penguin_tracks/processed/penguin_tracks_gemma_RAW_meta.csv'
# process
process.legacy_tracks(input_dir, output_dir, meta_out)
