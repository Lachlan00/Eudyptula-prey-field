# generate boat logs
source('scripts/CTD_functions.R')
setwd("/Users/lachlanphillips/Development/PhD/repos/Eudyptula")

# load the data
CTD_data <- cast.reader()
CTD_df <- CTD_data[[1]]
CTD_meta <- CTD_data[[2]]

#
ggplot(data=CTD_df, aes(x=temperature, y=depth, color=survey_id)) +
  geom_point(size=0.5) + scale_y_continuous(trans = "reverse") +
  facet_wrap(~survey_id)
  # theme(legend.position = "none")
  