## Get US census
#setwd('~/Dropbox/R/2017_seasonal_flu/')
source('~/Dropbox/R/ImmuneAgeStructure/census_sort.R')

# Import data
US.raw = read.csv('~/Dropbox/R/ImmuneAgeStructure/census_data_USA_1993-2015.csv', skip = 1, header = T)
years = unique(US.raw$Year)
# ----- Sort by year
US.demography = matrix(NA, nrow = length(years), ncol = 98, dimnames = list(years, rev(1918:2015)))
for(ii in 1:length(years)){
  year.pop = single.year.sort(US.raw[US.raw$Year == years[ii],])
  #For years < 2015, fill 2015:actual year with 0s and remove the same number of entries from end
  skips = 2015-years[ii]
  US.demography[ii, ] = c(rep(0, skips), year.pop[1:(98-skips)])
}
