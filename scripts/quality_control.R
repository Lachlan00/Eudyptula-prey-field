# Quality control functions
# Functions to assess data quality control and to
# filter data that does not meet set standards
library(ggplot2)

# List fo functions to add
# Dive detection
# points per track

#######################
# GPS fixes per track #
#######################
# Count number of non na lat values (so regularised tracks are counted)
quality.GPS.fix.count <- function(tracks, min_threshold=50, filter=F, plotMap=F){
  # split tracks by id
  tracks_split <- split(tracks, tracks$id)
  # for each track count number of fixes
  fixes <- sapply(tracks_split, function(x) sum(!is.na(x$lat)))
  # report
  if (plotMap){
    print(tracks.map(tracks, names(fixes)[which(fixes < min_threshold)],
                     title=paste0('GPS fix threshold < ',min_threshold)),legend=F)
  }
  plot(fixes, main='Total GPS fixes per track')
  abline(h=min_threshold, col='red', lty=2)
  fixes_low <- fixes
  fixes_low[fixes > min_threshold] <- NA
  points(fixes_low, col='red', add=T, pch=16)
  # FILTER
  if (filter){
    message('Filter is TRUE, dropping ',sum(!is.na(fixes_low)),' tracks..')
    tracks <- tracks[!tracks$id  %in% names(tracks_split)[which(fixes < min_threshold)],]
    tracks <- rows.and.levels(tracks)
    return(tracks)
  }
}

######################
# Gap to track ratio #
######################
# lat fixes that are NA (after gap filtering) vs having values
quality.gap.ratio <- function(tracks, max_threshold=1.2, filter=F, plotMap=F, specialPlot=F){
  # split tracks by id
  tracks_split <- split(tracks, tracks$id)
  # for each track computer the ratio of fixes to NAs
  gapRatio <- sapply(tracks_split, function(x) sum(is.na(x$lat))/sum(!is.na(x$lat)))
  # report
  if (plotMap){
    print(tracks.map(tracks, names(gapRatio)[which(gapRatio > max_threshold)],
                     title=paste0('Gap:Track Ratio Threshold > ',max_threshold)), legend=F)
  }
  plot(gapRatio, main='Gap:Track Ratio')
  abline(h=max_threshold, col='red', lty=2)
  gapRatio_high <- gapRatio
  gapRatio_high[gapRatio < max_threshold] <- NA
  points(gapRatio_high, col='red', add=T, pch=16)
  
  if (specialPlot){
    df <- as.data.frame(gapRatio)
    df$year <- sapply(tracks_split, function(x) year(x$dtUTC[1]))
    ggplot(df, aes(x=1:nrow(df), y=gapRatio, col=as.factor(year))) +
      geom_point() +
      labs(col='Year', x='index')
  }

  # FILTER
  if (filter){
    message('Filter is TRUE, dropping ',sum(!is.na(gapRatio_high)),' tracks..')
    tracks <- tracks[!tracks$id  %in% names(tracks_split)[which(gapRatio > max_threshold)],]
    tracks <- rows.and.levels(tracks)
    return(tracks)
  }
}

####################
# Incomplete trips #
####################
# Removes tracks that do not start and/or finish at the island
# threshold is in km
quality.incomplete <- function(tracks, start=T, end=T, filter=F, 
                               threshold=5, window=15, plotMap=F,
                               initCoords=c(150.2269, -36.25201)){
  # split tracks by id
  tracks_split <- split(tracks, tracks$id)
  # vector to hold results
  results <- c()
  # iterate through tracks and check
  for (i in 1:length(tracks_split)){
    track <- tracks_split[[i]]
    track <- track.displacement(track, initCoords=initCoords)
    results[i] <- sum(head(track, window)$disp < threshold*1000, na.rm=T) > 1 & 
                  sum(tail(track, window)$disp < threshold*1000, na.rm=T) > 1
  }
  # report
  filter.ls <- names(tracks_split)[!results]
  message(length(filter.ls),' tracks of ',length(tracks_split),' are incomplete')
  if (length(filter.ls)/length(tracks_split) > .5){
    warning('This is a lot of tracks to filter. This is likely due to tracks being ',
            'trimmed to first and last points on the ocean. Consider using a different ',
            'data source.')
  }
  # FILTER
  if (filter){
    message('Filter is TRUE, dropping ',length(filter.ls),' tracks..')
    tracks <- tracks[!tracks$id  %in% filter.ls,]
    tracks <- rows.and.levels(tracks)
    return(tracks)
  }
}

######################
# Filter by latitude #
######################
# I use to filter penguins going real far south (-36.5 min)
quality.coord.filter <- function(tracks, minlat=-90, maxlat=90, minlon=0, maxlon=360,
                                 plotMap=FALSE){
  # find min and max coordinates for penguin tracks
  min_coord <- aggregate(tracks[,c('lon','lat')], by=list(tracks$id), function(x) min(x, na.rm=T))
  names(min_coord) <- c('id', 'lon_min', 'lat_min')
  max_coord <- aggregate(tracks[,c('lon','lat')], by=list(tracks$id), function(x) max(x, na.rm=T))
  names(max_coord) <- c('id', 'lon_max', 'lat_max')
  minmax <- cbind(min_coord, max_coord[,c('lon_max', 'lat_max')])
  # filter by thresholds
  keep.ls <- minmax[minmax$lon_min > minlon & minmax$lat_min > minlat &
                  minmax$lon_max < maxlon & minmax$lat_max < maxlat,]
  # check if any filtering
  if (nrow(minmax) == nrow(keep.ls)){
    message('All tracks within bounds. No filtering was done.')
    return(tracks)
  } else {
    # removal list
    rm.ids <- minmax$id[!minmax$id %in% keep.ls$id]
    if (plotMap){
      print(tracks.map(tracks, as.character(rm.ids), legend=F))
    }
    # else make a subset
    track_sub <- tracks[tracks$id %in% keep.ls$id,]
    tracks_sub <- rows.and.levels(tracks_sub)
    return(tracks_sub)
  }
}
