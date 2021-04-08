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
model_dir <- './kriging/output/agg/models/'
plots_out_dir <- './kriging/output/agg/plots/'
quant_plots = FALSE
# set plot colours
cmap <- rev(cmocean('deep')(100))
clim_quanProb <- .95


# read in data
file_ls <- list.files(model_dir)
# filter for data bases only
file_ls <- file_ls[substr(file_ls, 11,12) == 'db']

# get variable list
var_ls <- unlist(lapply(file_ls, function(x) strsplit(x, '_')[[1]][3]))
vars <- "Sv_mean" # temp fix - need to be carefukl with those bloody underscores!

# cycle through variables
for (var in vars){
  message(paste0('\nProcessing ', var,'..'))
  # set long var name
  var_long <- paste0('Kriging.',var,'.estim')
  # get subset of the file list that is the variable
  file_idx <- which(var_ls == 'Sv') # NOTE!!! Another temp hack
  # load all the models for the variable
  message(paste('Loading',var,'models..'))
  df_ls <- pblapply(paste0(model_dir,file_ls[file_idx]), function(x) readRDS(x)@items)
  # get survey names
  survey_names <- lapply(file_ls[file_idx], function(s) strsplit(s,'_')[[1]][5:6])
  survey_names <- unlist(lapply(survey_names, function(s) substr(paste(s[1],s[2],sep='_'),1,7)))
  # find min and max value of variable across all data frames
  var_min <- min(unlist(lapply(df_ls, function(df) quantile(df[df$Polygon, var_long], probs=c(1-clim_quanProb), na.rm=T))))
  var_max <- max(unlist(lapply(df_ls, function(df) quantile(df[df$Polygon, var_long], probs=c(clim_quanProb), na.rm=T))))
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
    
    # TEMP hack! Override!
    # var_min <- min(df[df$Polygon, var_long],na.rm=T)
    # var_max <- max(df[df$Polygon, var_long],na.rm=T)
    
    # plotting
    #Set up dimension
    x <- unique(df$x1)
    y <- unique(df$x2)
    z <- -unique(df$x3)
    
    #Create Mesh
    #M <- mesh(x, y, z)
    #Create 3D data array
    p <- array(df[,var_long], dim=c(length(x),length(y),length(z)))
    # Hack P
    p[p < var_min] <- var_min + 1e-10
    p[p > var_max] <- var_max - 1e-10
    # Make closed 3D plot
    png(fn_3d_close, width=1000, height=850)
    slice3D(x/1000, y/1000, z, colvar = p,
            col = cmap,
            clim=c(var_min, var_max),
            NAcol="transparent",
            xs=NULL,
            ys = (y/1000)[c(1)],
            zs = c(z[c(1,21,41,61)]),
            d = 1,
            xlab="UTM Easting (km)",
            ylab="UTM Northing (km)",
            zlab="Depth (m)",
            clab=var,
            ticktype="detailed",
            main=paste('Survey', survey))
    dev.off()
    
    # hack data
    df[df[,var_long] < var_min & !is.na(df[,var_long]), var_long] <- var_min + 1e-10
    df[df[,var_long] > var_max & !is.na(df[,var_long]), var_long] <- var_max - 1e-10
    
    # cross section image
    xdepths <- c(0,5,10,15,20,25,30,35,40)
    xsection <- ggplot(data=df[df$x3 %in% xdepths,],
                       aes(x=x1/1000, y=x2/1000,fill=df[df$x3 %in% xdepths,var_long])) +
      geom_raster() +
      scale_fill_gradientn(colours=cmap, 
                           limits=c(var_min, var_max),na.value="#A0A0A0") +
      labs(fill=var) +
      ggtitle(paste('Survey', survey)) +
      xlab("UTM Easting (km)") +
      ylab("UTM Northing (km)") +
      facet_wrap(~x3) +
      theme(legend.key.height = unit(2, "cm"))
    # save plot
    ggsave(fn_xsection, xsection, width=9, height=8, dpi=200)
    # make qunatile plots
    if (quant_plots){
      fn_3d_quant_close <- paste0(plots_out_dir,"3dKrig_",var,"_quant_",survey,".png")
      fn_quant_xsection <- paste0(plots_out_dir,"3Dkrig_xsection_",var,"_quant_",survey,".png")
      #df[!df$Polygon, 'quants'] <- NA
      p <- array(as.numeric(df[,'quants']), dim=c(length(x),length(y),length(z)))
      # Make closed 3D plot
      # png(fn_3d_quant_close, width=1000, height=850)
      # slice3D(x/1000, y/1000, z, colvar = p,
      #         col = cmap,
      #         NAcol="transparent",
      #         xs=NULL,
      #         ys = (y/1000)[c(1)],
      #         zs = c(z[c(1,21,41,61)]),
      #         d = 1,
      #         xlab="UTM Easting (km)",
      #         ylab="UTM Northing (km)",
      #         zlab="Depth (m)",
      #         clab=paste(var,'quantiles'),
      #         ticktype="detailed",
      #         main=paste('Survey', survey))
      # dev.off()
      xsection <- ggplot(data=df[df$x3 %in% xdepths,],
                         aes(x=x1/1000, y=x2/1000,
                             fill=df[df$x3 %in% xdepths,'quants'])) +
        geom_raster() +
        scale_fill_manual(values=rev(cmocean('deep')(5))[2:5], na.translate = F) +
        labs(fill='Quantile') +
        ggtitle(paste('Survey', survey)) +
        xlab("UTM Easting (km)") +
        ylab("UTM Northing (km)") +
        facet_wrap(~paste(x3,'m')) +
        theme(legend.key.height = unit(2, "cm"))
      ggsave(fn_quant_xsection, xsection, width=9, height=8, dpi=200)
    }
  }
}

message('Complete!')
beep()
