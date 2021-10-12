require(tidyverse)

## Lach's data
d <- readRDS("data/pengData.rds")

## split into M and F datasets (simplest approach)
## recalculate diffs from mean weight
dm <- d %>%
  filter(sex == "m") %>%
  mutate(diffs1 = weight1 - mean(weight1)) %>%
  mutate(sex = "Male")

df <- d %>%
  filter(sex == "f") %>%
  mutate(diffs1 = weight1 - mean(weight1)) %>%
  mutate(sex = "Female")


## fit separate glm's for the 2 sexes
## remove intercept to test for differences from 0 (not grand mean)
m.glm <- glm(diffs1 ~ year - 1, data = dm)
summary(m.glm)
f.glm <- glm(diffs1 ~ year - 1, data = df)
summary(f.glm)

## plot data
d1 <- bind_rows(dm, df)

# make labels
# N labels
N.peng <- table(year(d$date))
N.peng <- data.frame(x=as.factor(2015:2019), 
                     y=278, #c(0,30,-20,42,-1), 
                     label=as.vector(N.peng),
                     m=as.vector(table(year(d$date[d$sex == 'm']))),
                     f=as.vector(table(year(d$date[d$sex == 'f']))))
N.peng$label.sex <- paste0('N=',N.peng$label,'\n(F',N.peng$f,':M',N.peng$m,')')

# stat labels
stat.lab <- data.frame(x=c(1:5 - .18, 1:5 + .18),
                       y=220,
                       sex=c(rep('Female',5), rep('Male',5)),
                       p=round(c(summary(f.glm)[["coefficients"]][,4],
                                 summary(m.glm)[["coefficients"]][,4]), 2))

p1 <- ggplot(d1, aes(year, diffs1, fill=sex)) +
  geom_text(data=N.peng,
            mapping=aes(x=x, y=y, label=label.sex), 
            inherit.aes = F) +
  geom_label(data=stat.lab,
             mapping=aes(x=x, y=y, label=p, fill=sex), 
             inherit.aes = F, alpha=.5, size=3) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_boxplot() +
  #facet_wrap(~ sex, ncol = 2) +
  theme(legend.position = "top") +
  ylab("Mass anomaly (g)") +
  xlab(element_blank()) +
  labs(fill='Sex') +
  theme_bw() +
  theme(text = element_text(size=18))

# get 2017 female % deviation from mean
grand.fmean <- mean(d1[d1$sex == 'Female','weight1'])
abs(mean(d1[d1$year == 2017 & d1$sex == 'Female','diffs1']))/grand.fmean*100

