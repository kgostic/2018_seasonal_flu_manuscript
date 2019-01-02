## This script imports raw data from AZ surveillance
## Script formats data for input into multinomial likelihood
## Also loads relevant probabilities of imprinting, and formats thses for input into likelihood
## Also generates age-specific indicators, which we use inside the likelihood to fit age specific step functions

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



names(raw.dat) = c('season', 'birthyear', 'subtype')
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

## Plots to check that agemats are properly formatted.
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

