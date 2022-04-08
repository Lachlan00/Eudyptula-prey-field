# Making 
#Calculate kernels
kde.df <- KDE_percentage(tracks[complete.cases(tracks),c('lon','lat')], 
                         levels=c(75, 95), facet=tracks[complete.cases(tracks),c('survey')])
kde.df$survey <- kde.df$facet
# Plots with contour lines

p.maps <- ggplot()  +  
  geom_raster(data = krig.dfs, aes(x=lon, y=lat, fill=Kriging.Sv_mean.estim)) +
  scale_fill_gradientn(colours=cmap, limits=c(vmin, vmax), na.value="grey66",
                       guide=guide_colorbar(barwidth=25, frame.colour=c("black"))) +
  theme_bw() + #grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL, fill=NULL) +
  theme(text = element_text(size=20)) +
  theme(legend.position = "bottom") +
  geom_path(data=kde.df, mapping=aes(x=x, y=y, group=group, linetype=factor(percent)),
            col='white') +
  geom_polygon(data=coast, mapping=aes(x=long, y=lat, group=group),
               fill='grey', col='#636363', size=.6, alpha=1, inherit.aes = F)+
  #geom_point(data=tracks, mapping=aes(x=lon, y=lat), col='darkgrey', alpha=.1) +
  facet_wrap(~survey) +
  coord_cartesian(xlim=xlim,  ylim=ylim, expand=c(0,0)) +
  guides(linetype = "none")

ggsave(p.maps, filename='output/paper_figures/simulation_paper/krigMap_pengs.png',
       bg="transparent",dpi=300, width=3500/300, height=3000/300)
