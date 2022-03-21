# Making 


# Plots with contour lines
p.maps <- ggplot()  +  
  geom_raster(data = krig.dfs, aes(x=lon, y=lat, fill=Kriging.Sv_mean.estim)) +
  scale_fill_gradientn(colours=cmap, limits=c(vmin, vmax), na.value="grey66",
                       guide=guide_colorbar(barwidth=25, frame.colour=c("black"))) +
  theme_bw() + #grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL, fill=NULL) +
  theme(text = element_text(size=18)) +
  theme(legend.position = "bottom") +
  facet_wrap(~survey) +
  stat_density_2d(data=tracks, mapping=aes(x=lon, y=lat), col='white', alpha=.5,
                  inherit.aes = F, size=.4, bins=10) +
  geom_polygon(data=coast, mapping=aes(x=long, y=lat, group=group),
               fill='grey', col='#636363', size=.6, alpha=1, inherit.aes = F)+
  coord_cartesian(xlim=xlim,  ylim=ylim, expand=c(0,0))

kernel_data <- ggplot_build(p.maps)$data[[2]]
kernel_data$survey <- c(2015:2019)[kernel_data$PANEL]
# find what levels are for each survey
kernel_data.ls <- split(kernel_data, kernel_data$survey)
for (i in 1:length(kernel_data.ls)){
  message(c(2015:2019)[i])
  print(unique(kernel_data.ls[[i]]$nlevel))
}
# manually adjust
kernel_data.ls$`2015` <- kernel_data.ls$`2015`[kernel_data.ls$`2015`$nlevel %in% c(0.25, 0.50, 0.75),]
kernel_data.ls$`2016` <- kernel_data.ls$`2016`[kernel_data.ls$`2016`$nlevel %in% c(0.25, 0.50, 0.75),]
kernel_data.ls$`2017` <- kernel_data.ls$`2017`[kernel_data.ls$`2017`$nlevel %in% c(0.2, 0.4, 0.8),]
kernel_data.ls$`2018` <- kernel_data.ls$`2018`[kernel_data.ls$`2018`$nlevel %in% unique(kernel_data.ls$`2018`$nlevel)[c(2,4,7)],]
kernel_data.ls$`2019` <- kernel_data.ls$`2019`[kernel_data.ls$`2019`$nlevel %in% c(0.2, 0.4, 0.8),]
kernel_data <- do.call(rbind, kernel_data.ls)

p.maps <- ggplot()  +  
  geom_raster(data = krig.dfs, aes(x=lon, y=lat, fill=Kriging.Sv_mean.estim)) +
  scale_fill_gradientn(colours=cmap, limits=c(vmin, vmax), na.value="grey66",
                       guide=guide_colorbar(barwidth=25, frame.colour=c("black"))) +
  theme_bw() + #grids(linetype = "dashed") +
  labs(x=NULL, y=NULL, alpha=NULL, size=NULL, fill=NULL) +
  theme(text = element_text(size=20)) +
  theme(legend.position = "bottom") +
  geom_path(data=kernel_data, mapping=aes(x=x, y=y, group=group), col='white',
                  inherit.aes = F, size=.4) +
  geom_polygon(data=coast, mapping=aes(x=long, y=lat, group=group),
               fill='grey', col='#636363', size=.6, alpha=1, inherit.aes = F)+
  facet_wrap(~survey) +
  coord_cartesian(xlim=xlim,  ylim=ylim, expand=c(0,0))

ggsave(p.maps, filename='output/paper_figures/simulation_paper/krigMap_pengs.png',
       bg="transparent",dpi=300, width=3500/300, height=3000/300)
