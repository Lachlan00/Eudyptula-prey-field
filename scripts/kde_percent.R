library(ggplot2)
library(ks)

# function to return contoour levels
KDE_percentage <- function(xy, levels=c(25, 50, 75), facet=NULL, plot=F){
  # if factet split data
  xy.orig <- xy
  if (!is.null(facet)){
    xy <- split(xy, facet)
    facets <- names(xy)
    xy.orig$facet <- facet
  } else {
    l <- list()
    l[[1]] <- xy
    xy <- l
  }
  
  # return the paths
  for (i in 1:length(xy)){
    kd <- kde(xy[[i]], compute.cont=TRUE)
    contours <- list()
    for (j in 1:length(levels)){
      contours[[j]] <- contourLines(x=kd$eval.points[[1]], y=kd$eval.points[[2]],
                                    z=kd$estimate, levels=kd$cont[paste0(100-levels[j], "%")])
      # sort gruoops for multi contours
      # if not x it is a list of list
      if (is.null(contours[[j]]$x)){
        for (k in 1:length(contours[[j]])){
          contours[[j]][[k]]$group <- k
          contours[[j]][[k]]$percent <- levels[j]
        }
        contours[[j]] <- lapply(contours[[j]], function(con) 
          data.frame(x=con$x, y=con$y, percent=con$percent,
                     level=con$level, group=con$group))
        contours[[j]] <- do.call(rbind, contours[[j]])
      } else {
        contours[[j]]$group <- 1
        contours[[j]]$percent <- levels[j]
        contours[[j]] <- data.frame(x= contours[[j]]$x, y= contours[[j]]$y, percent= contours[[j]]$percent,
                                    level= contours[[j]]$level, group= contours[[j]]$group)
      }
    }
    contours <- do.call(rbind, contours)
    xy[[i]] <- contours
    if (!is.null(facet)){
      xy[[i]]$facet <- facets[i]
    }
  }
  
  contour.df <- do.call(rbind, xy)
  contour.df$group <- paste0(contour.df$percent,'_',contour.df$group) 
  
  # plot to check
  if (plot){
    p <- ggplot() +
      geom_path(data=contour.df, mapping=aes(x=x, y=y, 
                                             linetype=factor(percent), 
                                             group=factor(group))) +
      geom_point(data=xy.orig, mapping=aes(x=xy.orig[,1], y=xy.orig[,2])) +
      labs(linetype="KDE%")
    if (!is.null(facet)){
      p <- p + facet_wrap(~facet)
    }
    print(p)
  }
  
  row.names(contour.df) <-NULL
  
  return(contour.df)
}

# Testing
# set.seed(420)
# N = 100
# df <- data.frame(lon=rnorm(N, mean=50, sd=25),
#                  lat = c(rnorm(N/2, mean=40, sd=15),
#                        rnorm(N/2, mean=60, sd=15)),
#                  state = rep(c('a','b'), each=N/2))
# 
# KDE_percentage(df[,c('lon','lat')], levels=c(75, 50, 90), plot=T)
# KDE_percentage(df[,c('lon','lat')], levels=c(10, 30, 50, 70), facet=df$state, plot=T)
              