# plot the 3D kriging fields

# libraries
library(pbapply)
library(plot3D)
library(plot3Drgl)
library(ggplot2)
library(cmocean)
library(beepr)
library(tools)

# setup
setwd("/Users/lachlanphillips/Development/PhD/repos/Eudyptula")
source('./scripts/visual.R')
model_dir <- './kriging/output/CTD/models/'
plots_out_dir <- './kriging/output/CTD/plots/'
# set plot colours
colour_maps <- list(cmaps=list(temperature = cmocean('thermal')(100),
                               salinity = cmocean('haline')(100),
                               soundVelocity = cmocean('speed')(100),
                               density = cmocean('dense')(100)),
                    label=list(temperature = 'Temperature (°C)',
                               salinity = 'Salinity (PSU)',
                               soundVelocity = 'Velocity (m/s)',
                               density = 'Density (kg/m^3)'))


# read in data
file_ls <- list.files(model_dir)
# filter for data bases only
file_ls <- file_ls[substr(file_ls, 11,12) == 'db']

# get variable list
var_ls <- unlist(lapply(file_ls, function(x) strsplit(x, '_')[[1]][3]))
vars <- unique(var_ls)

# cycle through variables
for (var in vars){
  message(paste0('\nProcessing ', var,'..'))
  # set long var name
  var_long <- paste0('Kriging.',var,'.estim')
  ############
  # TEMP FIX # Untill models rerun
  ############
  if (var_long == 'Kriging.soundVelocity.estim'){
    var_long <- 'Kriging.sound_velocity.estim'
  }
  ################
  # TEMP FIX END #
  ################
  # get subset of the file list that is the variable
  file_idx <- which(var_ls == var)
  # load all the models for the variable
  message(paste('Loading',var,'models..'))
  df_ls <- pblapply(paste0(model_dir,file_ls[file_idx]), function(x) readRDS(x)@items)
  # get survey names
  survey_names <- lapply(file_ls[file_idx], function(s) strsplit(s,'_')[[1]][4:5])
  survey_names <- unlist(lapply(survey_names, function(s) substr(paste(s[1],s[2],sep='_'),1,7)))
  
  # find min and max value of variable across all data frames
  var_min <- min(unlist(lapply(df_ls, function(df) min(df[df$Polygon, var_long],na.rm=T))))
  var_max <- max(unlist(lapply(df_ls, function(df) max(df[df$Polygon, var_long],na.rm=T))))
  # cycle thorugh and produce the plots
  message('Constructing plots..')
  for (i in 1:length(df_ls)){ 
    # load dataframe
    df <- df_ls[[i]]
    # outside Polygon gets NA
    df[!df$Polygon,var_long] <- NA
    # get survey name
    survey <- survey_names[i]
    # set file names
    # fn_3d_open <- paste0(plots_out_dir,"3dKrig_",var,"_",df$survey[1],"_open.png")
    fn_3d_close <- paste0(plots_out_dir,"3dKrig_",var,"_",survey,".png")
    fn_xsection <- paste0(plots_out_dir,"3Dkrig_xsection_",var,"_",survey,".png")
    
    # plotting
    #Set up dimension
    x <- unique(df$x1)
    y <- unique(df$x2)
    z <- -unique(df$x3)
    #Create Mesh
    #M <- mesh(x, y, z)
    #Create 3D data array
    p <- array(df[,var_long], dim=c(length(x),length(y),length(z)))
    # Make closed 3D plot
    png(fn_3d_close, width=1000, height=850)
    slice3D(x/1000, y/1000, z, colvar = p,
            col = colour_maps[['cmaps']][[var]],
            clim=c(var_min, var_max),
            NAcol="transparent",
            xs=NULL,
            ys = (y/1000)[c(1)],
            zs = c(z[c(1,21,41,61)]),
            d = 1,
            xlab="UTM Easting (km)",
            ylab="UTM Northing (km)",
            zlab="Depth (m)",
            clab=colour_maps[['label']][[var]],
            ticktype="detailed",
            main=paste('Survey', survey))
    dev.off()
    
    # cross section image
    xdepths <- c(0,10,20,30,50,80)
    xsection <- ggplot(data=df[df$x3 %in% xdepths,],
                       aes(x=x1/1000, y=x2/1000,fill=df[df$x3 %in% xdepths,var_long])) +
      geom_raster() +
      scale_fill_gradientn(colours=colour_maps[['cmaps']][[var]], 
                           limits=c(var_min, var_max),na.value="#A0A0A0") +
      labs(fill=colour_maps[['label']][[var]]) +
      ggtitle(paste('Survey', survey)) +
      xlab("UTM Easting (km)") +
      ylab("UTM Northing (km)") +
      facet_wrap(~paste(x3,'m')) +
      theme(legend.key.height = unit(2, "cm"))
    # save plot
    ggsave(fn_xsection, xsection, width=9, height=8, dpi=200)
  }
}

message('Complete!')
beep()
