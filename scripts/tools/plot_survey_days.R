# make transect plots for every day on the water
source('./scripts/CTD_functions.R')
source('./scripts/visual.R')

# __SETUP__
output_dir <- './output/survey_day_plots/'

# load all data
CTD_meta <- cast.reader()[[2]]

# split by survey
CTD_meta_ls <- split(CTD_meta, CTD_meta$survey_id)

# loop and plot
for (meta in CTD_meta_ls){
  # make the plot
  p <- plot.CTD.survey.path(meta)
  p <- p + ggtitle(meta$survey_id[1])
  ggsave(paste0(output_dir,meta$survey_id[1],'_surveydays.png') ,p)
}