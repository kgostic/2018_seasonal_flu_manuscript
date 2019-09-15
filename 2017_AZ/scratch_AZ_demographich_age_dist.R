## Check demographic age distribution over time in AZ

## setup workspace
setwd('~/Dropbox/R/2018_seasonal_flu/2017_AZ')
library('tidyverse')
library('ggplot2')
library('reshape')

## Load data
dat00_10 = read.csv('raw-data/census_by_state_2000_2010.csv')
dat10_18 = read.csv('raw-data/census_by_state_2010_2018.csv')
## Data from https://www2.census.gov/programs-surveys/popest/datasets/


# Extract AZ data, extract popestiamte columns
formatted_00_09 = dat00_10 %>% subset(NAME == 'Arizona' & SEX == 0 & AGE < 85) %>% select(c(NAME, AGE, contains('POPESTIMATE'))) %>% melt(id.vars = c('NAME', 'AGE')) %>% extract(col = 'variable', into = 'YEAR', regex = "\\w+(\\d\\d\\d\\d)") 
formatted_10_18 = dat10_18 %>% subset(NAME == 'Arizona' & SEX == 0 & AGE < 85) %>% select(c(NAME, AGE, contains('POPEST'))) %>% melt(id.vars = c('NAME', 'AGE')) %>% extract(col = 'variable', into = 'YEAR', regex = "POPEST(\\d\\d\\d\\d)_CIV") 

## Combine formatted data sets
pdat = rbind(subset(formatted_00_09, YEAR < 2010), 
             formatted_10_18)

## Add a column that counts the fraction in each age group
pdat = pdat %>% group_by(YEAR) %>% mutate(age_frac = value/sum(value))


## Plot
ggplot(pdat) +
  geom_line(aes(x = AGE, y = age_frac, color = YEAR))+
  ggtitle('Arizona Demography 2000 to 2018')+
  xlab('Age')+ylab('Population fraction')




