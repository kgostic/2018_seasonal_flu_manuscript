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
rm(list = ls(pattern = 'master'))
setwd('../2017_INSIGHT/')
source('00-Import_FLU002_-for-multinomial.R')

## Tidy INSIGHT
INSIGHTH1 = melt(rbind(H1.master), varnames = c('season', 'birth_year')); INSIGHTH1$subtype = 'H1'
INSIGHTH3 = melt(rbind(H3.master), varnames = c('season', 'birth_year')); INSIGHTH3$subtype = 'H3'
INSIGHTtidy = rbind(INSIGHTH1, INSIGHTH3); INSIGHTtidy$dataset = 'INSIGHT'
INSIGHTtidy=separate(INSIGHTtidy, col = season, into = c('country', 'season'), sep = '_')
head(INSIGHTtidy)


count.df$count[3:4] = c(sum(H1.master), sum(H3.master))


count.df


## Get total analyzed cases
sum(count.df$count)


## Get total number years analyzed, seasonal influenza
cs = c(rownames(H1.master), rownames(H3.master))
countries = unique(gsub(cs, pattern = '(\\w+)_\\w\\w.\\d?\\d?.?\\d\\d', replacement = '\\1')); countries
yrs_INSIGHT = 10:17


## All years, both data sets
all_yrs = unique(c(yrs_AZ, yrs_INSIGHT))

## Total years
length(all_yrs)

## Total countries
length(countries)



## Case counts from each data set, by subtype
write.csv(
AZtidy %>%
  group_by(subtype, season) %>%
  summarize(count = sum(value)), file = '../code_outputs/AZ_data_summary.csv', row.names = FALSE)

write.csv(
INSIGHTtidy %>%
  group_by(subtype, season) %>%
  summarize(count = sum(value)), file = '../code_outputs/INSIGHT_data_summary.csv', row.names = FALSE)

write.csv(
INSIGHTtidy %>%
  group_by(country) %>%
  summarize(count = sum(value)), file = '../code_outputs/INSIGHT_by_country.csv', row.names = FALSE)

