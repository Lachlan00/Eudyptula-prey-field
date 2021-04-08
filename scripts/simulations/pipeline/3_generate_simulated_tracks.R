# load functions
source('./scripts/simulations/simulate.R')

# load the data
model <- readRDS('./data/simulations/HMM_models/HMM.rds')
sim.df <- simulate.HMM(model, nbTracks=5000)

# save the data
saveRDS(sim.df, './data/simulations/sim_tracks.rds')

# make a pretty plot
tracks.map(sim.df, legend=F, plot_line=F, zoom=10)
