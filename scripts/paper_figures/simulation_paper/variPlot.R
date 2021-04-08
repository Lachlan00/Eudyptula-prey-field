# Load and plot variograms
# Load all
file.ls <- list.files('kriging/output/agg/models', pattern='vari*')
# keep only analysed 
file.ls <- file.ls[grepl('2015_S1|2016_S2|2017_S2|2018_S2|2019_S2',file.ls)]
# Load
vari.ls <- lapply(paste0('kriging/output/agg/models/',file.ls), function(x) readRDS(x))

# Plot
for (i in 1:length(vari.ls)){
  png(filename = paste0('kriging/output/agg/plots/vari_',(2015:2019)[i],'.png'),
      width = 600, height = 300)
  plot(vari.ls[[i]], main=paste((2015:2019)[i],'Variogram'))
  dev.off()
}

