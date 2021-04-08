library(ggplot2)
library(likert)
library(ggpubr)

# load data
analysis.ls <- readRDS('data/analysis_datasets/cumsum/simPaper_analysis.rds')
# transform into suitable dataset
boot.df <- do.call(c, lapply(analysis.ls, function(x) x$bootstrapped.sim.means))
boot.df <- data.frame(means=boot.df, year=factor(rep(2015:2019, each=100000)))
# get real means
real.means <- lapply(analysis.ls, function(x) x$consumption)
real.means <- lapply(real.means, function(x) mean(x[x$type == 'real','consumption']))
real.segments <- data.frame(x=seq(.5, 4.5, 1), xend=seq(1.5, 5.5, 1),
                            y=unlist(real.means), yend=unlist(real.means))

# Transform back to log
log.fun <- function(x){
  return(10*log10(x))
}

# Quick result report
stat.tests <- readRDS('./data/analysis_datasets/cumsum/simPaper_stats.rds')
print(data.frame(survey=unique(boot.df$year),
                 prop=unlist(lapply(stat.tests, 
                                    function(m) m$bootstrap.prop))))

# Make the plot
# p.boot <- ggplot(boot.df, aes(x=year, y=log.fun(means))) +
#   geom_violin(fill='grey85') +
#   geom_boxplot(width=0.1, fill='#2980b9', alpha=0.5) +
#   theme_bw() + grids(linetype = "dashed") +
#   theme(text = element_text(size=14)) +
#   theme(legend.position = "none") +
#   xlab(NULL) + ylab(bquote(~10^('Sv mean' / 10))) +
#   geom_segment(data=real.segments, mapping=aes(x=x, xend=xend, y=log.fun(y), yend=log.fun(yend)),
#                inherit.aes = F, col='red', linetype='dashed') +
#   theme(text = element_text(size=14))
#              
# ggsave(p.boot, filename='output/paper_figures/simulation_paper/bootPlot.png',
#        dpi=300, width=2500/300, height=2000/300)   
             
# calc diffs
boot.df.diff <- boot.df
boot.df.diff <- split(boot.df.diff, boot.df.diff$year)
for (i in 1:length(boot.df.diff)){
  boot.df.diff[[i]]$means <- log.fun(real.segments$y[i]) - log.fun(boot.df.diff[[i]]$means)
}
boot.df.diff <- do.call(rbind, boot.df.diff)
row.names(boot.df.diff) <- NULL

# Reverse  factor levels for coord flip (so 2015 is first to be plotted)
boot.df.diff$year <- reverse.levels(boot.df.diff$year)

# Attempt 2 (plot differences)
p.boot <- ggplot(boot.df.diff, aes(x=year, y=means)) +
  geom_violin(fill='#d35400', alpha=.6) +
  geom_boxplot(width=0.10, fill='lightgrey', alpha=0.93, outlier.shape=NA) +
  geom_hline(data=NULL, mapping=aes(yintercept=0), color='#2980b9', 
             linetype='dashed', size=1.08) +
  theme_bw() + grids(linetype = "dashed") +
  theme(text = element_text(size=14)) +
  theme(legend.position = "none") +
  xlab(NULL) + ylab('Search efficiency (Sv mean anomaly)') + 
  coord_flip()

ggsave(p.boot, filename='output/paper_figures/simulation_paper/bootPlot_pre.png',
       dpi=300, width=2500/300, height=2000/300)   
             
             