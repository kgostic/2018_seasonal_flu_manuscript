## Get total case counts from each data set
library(reshape2)
rm(list = ls())
setwd('~/Dropbox/R/2018_seasonal_flu/')
setwd('2017_AZ/')
source('00-Inputs_multinomial.R')
## Tidy AZ data
AZH1 = melt((H1.master), varnames = c('season', 'birth_year')); AZH1$subtype = 'H1'
AZH3 = melt((H3.master), varnames = c('season', 'birth_year')); AZH3$subtype = 'H3'
AZtidy = rbind(AZH1, AZH3); AZtidy$dataset = 'AZ'


yrs_AZ = c(93, 94, 95, 03, 04, 05, 06, 07, 10, 11, 12, 13, 14, 15)
count.df = data.frame(dataset = rep(c('AZ', 'INSIGHT'), each = 2),
                      subtype = rep(c('H1', 'H3'), times = 2),
                      kind = rep(c('seasonal', 'seasonal'), times = 2),
                      count = c(sum(H1.master), sum(H3.master), NA, NA))

## Case counts from each data set, by subtype
write.csv(
AZtidy %>%
  group_by(subtype, season) %>%
  summarize(count = sum(value)), file = '../tables/AZ_data_summary.csv', row.names = FALSE)