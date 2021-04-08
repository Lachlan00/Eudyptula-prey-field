# oad functions
source('./scripts/simulations/simulate.R')

# load the data
models <- sim.load.HMM.models()

# generate simulated penguin tracks
sim.tracks <- list()
for (i in 1:length(models)){
  message('\nSimulating tracks for ',names(models)[i])
  sim.tracks[[i]] <- simulate.HMM(models[[i]], nbTracks=50)
  sim.tracks[[i]]$year <- names(models)[i]
  sim.tracks[[i]]$survey_id <- models[[i]]$data$survey_id[1]
}
sim.df <- do.call(rbind, sim.tracks)

# save the data
saveRDS(sim.df, './data/simulations/sim_tracks.rds')

# make a pretty plot
tracks.map(sim.df, legend=F, colfactor='year', facet_year=T, plot_line=F)
