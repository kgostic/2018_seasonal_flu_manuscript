## This script imports raw data from AZ surveillance
## Script formats data for input into multinomial likelihood
## Also loads relevant probabilities of imprinting, and formats thses for input into likelihood
## Also generates age-specific indicators, which we use inside the likelihood to fit age specific step functions
library('dplyr')
library('reshape')

## Load Arizona data
raw.dat = read.csv('raw-data/AZ_flu_surveillance_93-94.csv')
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_94-95.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_02-03.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_03-04.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_04-05.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_05-06.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_06-07.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_07-08.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_08-09.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_09-10.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_10-11.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_11-12.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_12-13.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_13-14.csv'))
raw.dat = rbind(raw.dat, read.csv('raw-data/AZ_flu_surveillance_14-15.csv'))
## Exclude one case whose birth year is mis-entered
raw.dat = raw.dat[-364, ]
## Calculate age for plotting
raw.dat$year = as.numeric(gsub(raw.dat$SEASON, pattern = '(\\d{4})\\d{2}', replacement = '\\1'))
raw.dat$age = raw.dat$year+1 - raw.dat$BIRTHYEAR
names(raw.dat) = c('season', 'birthyear', 'subtype', 'year', 'age')
## Exclude cases born before 1918 (imprinting history not known, and there are few cases)
raw.dat = raw.dat[-which(raw.dat$birthyear < 1918), ]
write.csv(raw.dat, file = 'processed-data/AZ_seasonal_linelist.csv', row.names = FALSE)


## Tabluate the number of H1 and H3 cases in each birth year where cases were observed
ssns = unique(raw.dat$season)
bys = 2015:1918
H1.master = H3.master = matrix(NA, nrow = length(ssns), ncol = length(bys), dimnames = list(ssns, bys))
for(ss in ssns){
  validH1 = raw.dat[which(raw.dat$season == ss & raw.dat$subtype == 'H1'), ]
  validH3 = raw.dat[which(raw.dat$season == ss & raw.dat$subtype == 'H3'), ]
  for(ii in 1:length(bys)){
    H1.master[as.character(ss), ii] = sum(validH1$birthyear == bys[ii])
    H3.master[as.character(ss), ii] = sum(validH3$birthyear == bys[ii])
  }
}




## Reconstruct imprinting patterns
## Commented code was used to generate reconstructions
## Instead of re-running this code each time, we can just load the saved outputs
## Commented code is provided for reproducibility
# setwd('../../Reconstructions/')
# source('01-reconstruct-imprinting-histories.R')
# weights = get.weights.master(years.out = c(1994, 1995, 2003:2015), Countries.out = 'USA', region.in = c('default'))
# save(weights, file = paste('AZ_weights_', Sys.Date(), '.RData', sep = ''))
# setwd('../2018_seasonal_flu/2017_AZ/')
load('../../Reconstructions/AZ_weights_2018-11-29.RData')
## Extract relevant protection status from weights outputs.
## Loaded variable is called weights, and is a list of length 4
## weights[[1]] gives probabilities of childhood imptinting to H1N1
## weights[[2]] gives probabilities of imprinting to H2N2
## weights[[3]] gives probs of H3N2 imrinting
## weights[[4]] gives probabilities of remaining naive. This only takes non-zero probabilities in cohorts under age 13.
proH1.master = weights[[1]][,as.character(2015:1918)] # Individuals who imprinted to H1N1 will be protected against H1N1 at the HA subtype level
proH2.master = weights[[2]][,as.character(2015:1918)]
proH3.master = weights[[3]][,as.character(2015:1918)] # Individuals who imprinted to H3N2 will be protected against H3N2 at the HA subtype level
prog1.master = weights[[1]][,as.character(2015:1918)]+weights[[2]][,as.character(2015:1918)] # Individuals who imprinted to H1N1 OR H2N2 will be protected against H1N1 at the HA group level
prog2.master = proH3.master # Individuals who imprinted to H3N2 will be protected against H3N2 at the HA group level
proN1.master = proH1.master # Individuals who imprinted to H1N1 will be protected against H1N1 at the NA subtype level
proN2.master = weights[[2]][,as.character(2015:1918)]+weights[[3]][,as.character(2015:1918)] # Individuals who imprinted to H2N2 or H3N2 will be protected against H3N2 at the NA subtype level




#### Define age-specific indicators so you can convert between birthyear and age in different years of case observation
agemat = t(sapply(X = c(1994, 1995, 2003:2015), function(xx) xx - 2015:1918))
rownames(agemat) = ssns; colnames(agemat) = bys

## Define age bins to which we will fit a setp function
a0.4 = a5.10 = a11.17 = a18.24 = a25.31 = a32.38 = a39.45 = a46.52 = a53.59 = a60.66 = a67.73 = a74.80 = a81.90plus = H1.master*0 # Initialize
a0.4[agemat %in% 0:4] = 1
a5.10[agemat %in% 5:10] = 1
a11.17[agemat %in% 11:17] = 1
a18.24[agemat %in% 18:24] = 1
a25.31[agemat %in% (25:31)] = 1
a32.38[agemat %in% (32:38)] = 1
a39.45[agemat %in% (39:45)] = 1
a46.52[agemat %in% (46:52)] = 1
a53.59[agemat %in% (53:59)] = 1
a60.66[agemat %in% (60:66)] = 1
a67.73[agemat %in% (67:73)] = 1
a74.80[agemat %in% (74:80)] = 1
a81.90plus[agemat >= 81] = 1
## Each agemat is an indicator matrix with 1s to indicate which birth years belond to a certain age class in the year of case observation





## Demography matrix
## Load data
dat00_10 = read.csv('raw-data/census_by_state_2000_2010.csv')
dat10_18 = read.csv('raw-data/census_by_state_2010_2018.csv')
## Data from https://www2.census.gov/programs-surveys/popest/datasets/
# Extract AZ data, extract popestiamte columns
formatted_00_09 = dat00_10 %>% subset(NAME == 'Arizona' & SEX == 0 & AGE <= 85) %>% select(c(NAME, AGE, contains('POPESTIMATE'))) %>% melt(id.vars = c('NAME', 'AGE')) %>% tidyr::extract(col = variable, into = 'YEAR', regex = "\\w+(\\d\\d\\d\\d)") 
formatted_10_18 = dat10_18 %>% subset(NAME == 'Arizona' & SEX == 0 & AGE <= 85) %>% select(c(NAME, AGE, contains('POPEST'))) %>% melt(id.vars = c('NAME', 'AGE')) %>% tidyr::extract(col = 'variable', into = 'YEAR', regex = "POPEST(\\d\\d\\d\\d)_CIV") 

## Combine formatted data sets
pdat = rbind(subset(formatted_00_09, YEAR < 2010),
             formatted_10_18)

## Write a function to reformat into standard model input format
rys = as.numeric(gsub(rownames(H1.master), pattern = '(\\d{4}).+', replacement = "\\1")) # extract first year from rownames
demog = H1.master*0 # Set up empty matrix
## From demographic data pdf: POPESTIMATE2000 7/1/2000 resident population estimate. Because popest starts in July, use 2000 for the 2000-2001 season
for(rowindex in 1:nrow(demog)){
  rawdem = subset(pdat, YEAR == max(rys[rowindex], 2000))$value # gives a vector of # in age from 0 to 84, last entry gives total for 85+
  ## Convert from age to birth year using age mat
  startcol = which(agemat[rowindex,]==0)
  ## Age class "85" is really the total for all adults 85 and older. 
  ## We could just evenly distribute this total count over ages 85:100, but the data show that the number of people in each age class tends to decrease with increasing age. So even distrbution would systematically overestimate the number of 100-year-olds, and underestimate the number of 85-year-olds.
  ## Instead, fit a linear model to the shape of decrease from ages 75-84, the oldest 10 years in the single-year-of-age data
  pdf = data.frame(xx = 75:84, yy = rawdem[76:85]) ## Extract data to which we fit
  ft = lm(formula = yy~xx, data = pdf) ## Fit a linear model
  pd = predict(ft, data.frame(xx = c(85:97))) ## Use lm to extrapolate number of individuals in each age class
  pd[which(pd<1700)] = 1000 # Set minimum number in any single year of age.
  ## Check interpolation
  # plot(0:84, orig[1:85], xlim = c(0, 99), ylim = c(0, 90000))
  # points(85:97, pd, col = 'blue')
  # sum(pd); rawdem[86] # Should be similar in value
  rawdem[86:98] = pd ## Fill in interpolation
  #points(0:97, rawdem, pch = 16, cex = .8, col = 'red')
  demog[rowindex, startcol:98] = rawdem[1:sum(agemat[rowindex,]>=0)]
}
## normalize
demog = demog/rowSums(demog)
## Check inputs by birth year
# cols = rainbow(16)
# plot(2015:1918, demog[1,], col = cols[1])
# for(ii in 2:nrow(demog)){
#   lines(2015:1918, demog[ii,], col = cols[ii])
# }
# ## Plots to check that agemats are properly formatted.
# barplot(a18.24)
# barplot(a25.31)
# barplot(a32.38)
# barplot(a39.45)
# barplot(a46.52)
# barplot(a53.59)
# barplot(a60.66)
# barplot(a67.73)
# barplot(a74.80)
# barplot(a81.90plus)
#barplot(a0.4 + a5.10 + a11.17 + a18.24+a25.31+a32.38+a39.45+a46.52+a53.59+a60.66+a67.73+a74.80+a81.90plus, col = rainbow(14))


## Get a table of counts to input into table S1
table(raw.dat$season)


## Hold out the 2009 pandemic (2008-09 and 2009-10 seasons)
a0.4_2009 = a0.4[c('200809', '200910'), ]
a5.10_2009 = a5.10[c('200809', '200910'), ]
a11.17_2009 = a11.17[c('200809', '200910'), ]
a18.24_2009 = a18.24[c('200809', '200910'), ]
a25.31_2009 = a25.31[c('200809', '200910'), ]
a32.38_2009 = a32.38[c('200809', '200910'), ]
a39.45_2009 = a39.45[c('200809', '200910'), ]
a46.52_2009 = a46.52[c('200809', '200910'), ]
a53.59_2009 = a53.59[c('200809', '200910'), ]
a60.66_2009 = a60.66[c('200809', '200910'), ]
a67.73_2009 = a67.73[c('200809', '200910'), ]
a74.80_2009 = a74.80[c('200809', '200910'), ]
a81.90plus_2009 = a81.90plus[c('200809', '200910'), ]
H1.master_2009 = H1.master[c('200809', '200910'), ]
prog1.master_2009 = prog1.master[c('2009USA', '2010USA'), ]
proH1.master_2009 = proH1.master[c('2009USA', '2010USA'), ]
proN1.master_2009 = proN1.master[c('2009USA', '2010USA'), ]
H3.master_2009 = H3.master[c('200809', '200910'), ]
prog2.master_2009 = prog2.master[c('2009USA', '2010USA'), ]
proH3.master_2009 = proH3.master[c('2009USA', '2010USA'), ]
proN2.master_2009 = proN2.master[c('2009USA', '2010USA'), ]
demog_2009 = demog[c('200809', '200910'), ]


## Remove 2009-2010 season from master data for model fitting
a0.4 = a0.4[-c(9,10), ]
a5.10 = a5.10[-c(9,10), ]
a11.17 = a11.17[-c(9,10), ]
a18.24 = a18.24[-c(9,10), ]
a25.31 = a25.31[-c(9,10), ]
a32.38 = a32.38[-c(9,10), ]
a39.45 = a39.45[-c(9,10), ]
a46.52 = a46.52[-c(9,10), ]
a53.59 = a53.59[-c(9,10), ]
a60.66 = a60.66[-c(9,10), ]
a67.73 = a67.73[-c(9,10), ]
a74.80 = a74.80[-c(9,10), ]
a81.90plus = a81.90plus[-c(9,10), ]
H1.master = H1.master[-c(9,10), ]
prog1.master = prog1.master[-c(9,10), ]
proH1.master = proH1.master[-c(9,10), ]
proN1.master = proN1.master[-c(9,10), ]
H3.master = H3.master[-c(9,10), ]
prog2.master = prog2.master[-c(9,10), ]
proH3.master = proH3.master[-c(9,10), ]
proH2.master = proH2.master[-c(9,10), ]
proN2.master = proN2.master[-c(9,10), ]
demog = demog[-c(9,10), ]



