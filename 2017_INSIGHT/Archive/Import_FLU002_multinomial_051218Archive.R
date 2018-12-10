## This file loads FLU002 data for analysis.
setwd('~/Dropbox/R/2017_INSIGHT/')
rm(list = ls()) # Clear memory


# Load data
dat.002 = read.csv('data002.csv', colClasses = c('character', 'integer', 'integer', 'integer', 'character', 'integer', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character'))

# Enter actual country names instead of country codes
Country = character(nrow(dat.002))
ccs = unique(dat.002$COUNTRY_CODE); ccs # Existing codes
cnames = c('Denmark', 'Spain', 'Germany', 'Estonia', 'USA', 'Belgium', 'Portugal', 'Poland', 'Austria', 'UK', 'Australia', 'Thailand', 'Argentina', 'Chile', 'Greece', 'Peru', 'Japan') # Names with which to replace 
rbind(ccs, cnames) # Check
for(ii in 1:length(ccs)){
  Country[which(dat.002$COUNTRY_CODE == ccs[ii])] = cnames[ii]
}
dat.002$country = Country
rm(cnames, ccs, Country)

# Exclude cases over age 90
# Do this because individuals observed in 2009 of age 90 could have been born in 1918
# Older individuals could have been born before 1918
dat.002 = dat.002[-which(dat.002$age > 90), ]

# Make sure all month entries follow the specified format
if(any(!(dat.002$month %in% c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov','Dec')))){ error('Ivalid month passed to month_to_num')}
month_to_num = function(months.in){
  months.in[months.in == 'Jan'] = 1
  months.in[months.in == 'Feb'] = 2
  months.in[months.in == 'Mar'] = 3
  months.in[months.in == 'Apr'] = 4
  months.in[months.in == 'May'] = 5
  months.in[months.in == 'Jun'] = 6
  months.in[months.in == 'Jul'] = 7
  months.in[months.in == 'Aug'] = 8
  months.in[months.in == 'Sep'] = 9
  months.in[months.in == 'Oct'] = 10
  months.in[months.in == 'Nov'] = 11
  months.in[months.in == 'Dec'] = 12
  as.numeric(months.in)
}
# Convert
dat.002$num_month = month_to_num(dat.002$month)

# Calculate imprinting probs based on birth year
table(dat.002$COUNTRY_CODE, dat.002$year) # View countires and years in which we need to do reconstrutions
source('../Reconstructions/Infection_age_structure_2017.R')
# This code used to generate and save a list of weights. 
# Weights are saved and loaded in future runs to save runtime.
# weights = get.type.weights.AB(years.out = 2009:2017, Countries.out = c('Argentina', 'Australia', 'Austria', 'Belgium', 'Chile', 'Germany', 'Denmark', 'Spain', 'Estonia', 'UK', 'Greece', 'Japan', 'Peru', 'Poland', 'Portugal', 'Thailand', 'USA'), region.in = c('defulat', 'default', 'Euro', 'Euro', 'defualt', 'Euro', 'Euro', 'Euro', 'Euro', 'Euro', 'Euro', 'Asia', 'default', 'Euro', 'Euro', 'Asia', 'default'))
# save(weights, file = 'INSIGHT002_weights.RData')
load('INSIGHT002_weights.RData')

# Define seasons based on October to October years
sort(unique(dat.002$year))
NHfirst.part = c('Oct', 'Nov', 'Dec')
NHsecond.part = c('Jan', 'Feb', 'Mar')
SH = c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
ssn = vector('character', nrow(dat.002)) # Create a new vector to record season codes
season.codes.NH = paste('NH', c('09', 10:16), 10:17, sep = '.')
season.codes.SH = paste('SH', 10:17, sep = '.')
# Assign each observation to a NH or SH season (Oct-March = NH, April - Sept = SH)
yy = 2009
for(ii in 1:8){
  NH.indices = c(which(dat.002$year == yy & dat.002$month %in% NHfirst.part), which(dat.002$year == yy+1 &  dat.002$month %in% NHsecond.part))
  ssn[NH.indices] = season.codes.NH[ii]
  SH.indices = which(dat.002$year == yy + 1 & dat.002$month %in% SH)
  ssn[SH.indices] = season.codes.SH[ii]
  yy = yy + 1
}
# Add the season assignements to the data frame
dat.002$season = (ssn)
rm(ssn)
# 
full.dat = dat.002 # Archive all columns

#####################################################
## Pull out variables used in model
#####################################################
dat.002 = subset(dat.002, select = c('age', 'anyvac', 'anydx', 'anyav', 'season', 'country', 'season', 'flutype'))
which(apply(dat.002, MARGIN = 1,  function(xx) (any(is.na(xx)))))
# Exclude rows with any NAs
dat.002 = na.omit(dat.002)
# Exclude those with unknown vaccination status.
invalid = which(dat.002$anyvac == 2)
dat.002 = dat.002[-invalid, ]
# Exclude two observations with very low age
invalid = which(dat.002$age == 17)
dat.002 = dat.002[-invalid,]

# check case counts by country and subtype
table(dat.002$country, dat.002$flutype) -> tab1
tab1
tab1[,1:2]
rowSums(tab1[,1:2])

## Exclude cases with unknown flutype or coinfection
dat.002 = subset(dat.002, flutype %in% c(1,2,3,4))

# Null construction 1 (balanced number of cases from H1, H3, B, neg):
#  From each country, sample max(125, min_of_4_categories) from each category.
#  Use resulting age distribution to formulate null
## Analyze data from Argentina, Thailand and Belgium. These countries meet the following inclusion criteria:
# * At least 100 cases observed in the country in each of the following categories:
#     - H1N1 confirmed
#     - H3N2 confirmed
#     - FluB confirmed
#     - Flu negative
null.inds = NULL
ccs = unique(dat.002$country)
for(cc in ccs){ # for each country
  valid = subset(dat.002, country == cc)
  if(!any(table(valid$flutype) < 150)){ # If all flutypes: H1, H3, B and Neg all have >=100 observations in the country of interest, then include these indices.
    counts = table(valid$flutype) #
    min.count = min(counts) # Get the minimum number of observations in any bin
    for(ftype in 1:4){
      typeInds = which(dat.002$country == cc & dat.002$flutype == ftype) # Get the indices of all rows with the correct country and flutype
      sampInds = sample(x = typeInds, size = min.count, replace = FALSE) # Randomly sample a number equal to min.count of these indices
      null.inds = c(null.inds, sampInds) # Append the sampled indices to the overall vector of null distribution indices
    } # End loop over flutype 1:4
  } # End if statement (move on if not enough obs in the country of interest)
} # End loop over country
null_set_1 = dat.002[null.inds, ] # Extract desired rows from the full data
table(null_set_1$country, null_set_1$flutype) # Check that null set provides an even sample






# Null construction 2 (even H1, H3 + all others):
# Null construction 1 (balanced number of cases from H1, H3, B, neg):
#  From each country, sample max(100, min_of_4_categories) from each category.
#  Use resulting age distribution to formulate null
## Analyze data from Argentina, Thailand, USA and Belgium. These countries meet the following inclusion criteria:
# * At least 150 cases observed in the country in each of the following categories:
#     - H1N1 confirmed
#     - H3N2 confirmed
#     - FluB confirmed
#     - Flu negative
#  From each country, sample max(100, min_of_4_categories) from H1N1 and H3N2
#  Include all fluB and null cases
#  Use resulting age distribution to formulate null
null.inds = NULL
ccs = unique(dat.002$country)
for(cc in ccs){ # for each country
  valid = subset(dat.002, country == cc)
  if(!any(table(valid$flutype)[1:2] < 150)){ # If flutypes: H1, H3, have >=100 observations in the country of interest, and flub + flu neg have >200 then include these indices.
    counts = table(valid$flutype)[1:2] #
    min.count = min(counts) # Get the minimum number of observations in any bin
    for(ftype in 1:2){
      typeInds = which(dat.002$country == cc & dat.002$flutype == ftype) # Get the indices of all rows with the correct country and flutype
      sampInds = sample(x = typeInds, size = min.count, replace = FALSE) # Randomly sample a number equal to min.count of these indices
      null.inds = c(null.inds, sampInds) # Append the sampled indices to the overall vector of null distribution indices
    } # End loop over flutype 1:4
    # For flu B and neg, use all obs from the country of interest
    typeInds = which(dat.002$country == cc & dat.002$flutype %in% c(3,4)) # Get the indices of all rows with the correct country and flutype
    null.inds = c(null.inds, typeInds) # Append the sampled indices to the overall vector of null distribution indices
  } # End if statement (move on if not enough obs in the country of interest)
} # End loop over country
null_set_2 = dat.002[null.inds, ] # Extract desired rows from the full data
table(null_set_2$country, null_set_2$flutype) # Check that null set provides an even sample of H1, H3



 
# Null construction 3 recommended by Jamie (Use all B and neg cases):
# Inclusion criteria:
#    * At least 100 H1N1 and 100 H3N2 cases observed in the country of interest
#    * At least 400 cases of flu B or PCR negative observed in the country of interest
## Analyze data from Argentina, Thailand, USA and Belgium. These countries meet the following inclusion criteria:
# * At least 150 cases observed in the country in each of the following categories:
#     - H1N1 confirmed
#     - H3N2 confirmed
#     - FluB confirmed
#     - Flu negative
#  From each country, sample max(100, min_of_4_categories) from H1N1 and H3N2
#  Include all fluB and null cases
#  Use resulting age distribution to formulate null
null.inds = NULL
ccs = unique(dat.002$country)
for(cc in ccs){ # for each country
  valid = subset(dat.002, country == cc)
  tab = table(valid$flutype)
  if(tab[1]>=100&tab[2]>=100&(tab[3]+tab[4])>=400){ # If flutypes: H1, H3, have >=100 observations in the country of interest, and flub + flu neg have >200 then include these indices.
      newInds = which(dat.002$country == cc & dat.002$flutype %in% c(3,4)) # Get the indices of all rows with the correct country and flutype B or PCR Neg
      null.inds = c(null.inds, newInds) # Append the sampled indices to the overall vector of null distribution indices
  } # End if for inclusion of country
}# End loop over country
null_set = dat.002[null.inds, ] # Extract desired rows from the full data
table(null_set$country, null_set$flutype) # Check that null set provides an even sample of H1, H3







adH1 = subset(dat.002, flutype == 1, select = age)$age
adH1 = hist(adH1, breaks = seq(17.5, 90.5, by = 1))$counts
adH3 = subset(dat.002, flutype == 2, select = age)$age
adH3 = hist(adH3, breaks = seq(17.5, 90.5, by = 1))$counts

par(mfrow = c(1,1))
ad1 = hist(null_set_1$age, breaks = seq(17.5, 90.5, by = 1))$counts
ad2 = hist(null_set_2$age, breaks = seq(17.5, 90.5, by = 1))$counts
ad = hist(null_set$age, breaks = seq(17.5, 90.5, by = 1))$counts
#barplot(rbind(ad1, ad2)/c(sum(ad1), sum(ad2)), beside = T)
plot(18:90, ad1/sum(ad1), pch = 16, col = 'tan', cex = .7, ylab = 'fraction', ylim = c(0, .03), xlab = 'age')
lines(density(x = rep(18:90, ad1)), col = 'tan')
#lines(18:90, ad1/sum(ad1), col = 'tan')
points(18:90, ad2/sum(ad2), pch = 16, col = 'gray', cex = .7)
lines(density(x = rep(18:90, ad2)), col = 'gray')
# lines(18:90, ad2/sum(ad2), col = 'gray')
points(18:90, ad/sum(ad), pch = 16, col = 'black', cex = .7)
lines(density(x = rep(18:90, ad)), col = 'black')
# points(18:90, adH1/sum(adH1), pch = 16, col = 'dodgerblue', cex = .7)
# lines(density(x = rep(18:90, adH1)), col = 'dodgerblue')
# #lines(18:90, adH1/sum(adH1), col = 'dodgerblue')
# points(18:90, adH3/sum(adH3), pch = 16, col = 'red', cex = .7)
# #lines(18:90, adH3/sum(adH3), col = 'red')
# lines(density(x = rep(18:90, adH3)), col = 'red')
legend('topright', c('All balanced null', 'Balanced H1, H3 null', 'B and NEG only', 'Confirmed H1', 'Confirmed H3', 'smoothed density'), col = c('tan', 'gray', 'black', 'dodgerblue', 'red', 'black'), pch = c(16, 16, 16, 16, 16, NA), lty = c(NA, NA, NA, NA, NA, 1))



### USE DENSITY TO GET AGE NULL
### USE SOME OTHER SMOOTHER TO GET VAC, DX, ETC.
########################
## This is an exercise to verify that loess can safely be used to generate smoothed 
##    proportion estimates from binary data.
## This exercise also aims to verify that if input y values are a sample proportion estimate for the 
##    age of interest, that loess weights should be input as a vector: (n1, n2, n3, ...)
##    and not as (1/n1, 1/n2, 1/n3...)
## First, I define a set of weights (numbers of observations at each age)
{
ARGEN = subset(dat.002, country == "Argentina")
wts = table(ARGEN$age)
xshort = as.numeric(names(wts))
## Create an expanded x vector with each age value repeated the number of times specified by the weights vector
xlong = rep(xshort, wts)
## Sample binary y values corresponding to x values with probability of drawing a 1 = 0.01*x
## ylong stores a binary observation similar to the observed vaccination data. With 0s and 1s repeated for multiple observations at each age.
## yshort summarizes the sample proportion of 1s at each age value
probs = (xshort-17)*.001
ylong = NULL
yshort = numeric(length(xshort))
for(ii in 1:length(xshort)){
  ynew = sample(x = c(1,0), size = wts[ii], replace = TRUE, prob = c(probs[ii], 1-probs[ii]))
  ylong = c(ylong, ynew)
  yshort[ii] = mean(ynew)
}
## Now plot the long and short data sets.
plot(jitter(xlong), ylong)
points(xshort, yshort, col = 'red')
## Create data frames (makes plotting easier)
longdat = data.frame(x = xlong, y = ylong)
shortdat = data.frame(x = xshort, y = yshort)
## Fit loess smoothers to the long and short data
longfit = loess(formula = y ~ x, data = longdat, degree = 1)
shortfit = loess(formula = y~x, data = shortdat, weights = wts, degree = 1)
# Plot the fits to verify similar outcomes with either strategy
lines(xshort, predict(longfit, data.frame(x = xshort)))
lines(xshort, predict(shortfit, data.frame(x = xshort)), col = 'red')
lines(xshort, (xshort-17)*.001, col = 'blue')
legend(25, .8, c('Simulated binary data', 'Simulated summary data', 'loess fit to binary data', 'loess fit to summary data', 'true prob relationship'), col = c('black', 'red', 'black', 'red', 'blue'), pch = c(1,1,NA,NA,NA), lty = c(NA,NA,1,1,1), bty = 'n')
}
## Looks reasonable.
## I used the actual numbers of individuals of each age from Argentina, so these fits are representative of how well we might expect to approximate the true population pattern using this smoother.
## Most of the time, the soother approximates the true, underlying population relationship.
## The smoothed estimate doesn't always match the underlying population pattern perfectly, but it appears to be an unbiased estimator (I've run this several times to get a feel for the range of possible outcomes).
## The soothed estimate from the binary data usually matches the true undelrying pattern a bit better than the weighted summary smoother, so I'll use the raw, binary data to fit my model inputs
##### End exercise ##################
######################################



### USE DENSITY TO GET AGE NULL
### USE SOME OTHER SMOOTHER TO GET VAC, DX, ETC.
# ########################
# ## Repeat the above exercise 100 times and plot the outcomes
# pdf('cross_validation1.pdf')
# par(mfrow = c(2,2))
# tnsblack = rgb(red = 0,green = 0,blue = 0,alpha = 100, maxColorValue = 255)
#   ARGEN = subset(dat.002, country == "Argentina")
#   wts = table(ARGEN$age)
#   xshort = as.numeric(names(wts))
#   ## Create an expanded x vector with each age value repeated the number of times specified by the weights vector
#   xlong = rep(xshort, wts)
# 
# ## Use cross validation to choose span
#   spans = seq(.05, 1.5, by = .05)
#   ssse = numeric(length(spans))
#   for(ss in 1:length(spans)){
#   ## Sample binary y values corresponding to x values with probability of drawing a 1 = 0.01*x
#   ## ylong stores a binary observation similar to the observed vaccination data. With 0s and 1s repeated for multiple observations at each age.
#   ## yshort summarizes the sample proportion of 1s at each age value
#   probs = -.0001*(xshort-50)^2+.19 #(xshort-17)*.001
#   plot(xshort, probs, col = 'blue', ylim = c(0, 1), type = "l", main = paste('span=', spans[ss]))
#   preds = matrix(NA, nrow = 100, ncol = length(xshort))
#   sse = numeric(100)
#   for(rr in 1:100){
#   ylong = NULL
#   yshort = numeric(length(xshort))
#   for(ii in 1:length(xshort)){
#     ynew = sample(x = c(1,0), size = wts[ii], replace = TRUE, prob = c(probs[ii], 1-probs[ii]))
#     ylong = c(ylong, ynew)
#     yshort[ii] = mean(ynew)
#   }
#   ## Create data frames (makes plotting easier)
#   longdat = data.frame(x = xlong, y = ylong)
#   shortdat = data.frame(x = xshort, y = yshort)
#   ## Fit loess smoothers to the long and short data
#   longfit = loess(formula = y ~ x, data = longdat, degree = 2, span = spans[ss])
#   #shortfit = loess(formula = y~x, data = shortdat, weights = wts)
#   # Plot the fits to verify similar outcomes with either strategy
#   preds[rr, ] = predict(longfit, data.frame(x = xshort))
#  # lines(xshort, preds[rr, ])
#   se = sum((preds[rr,]-probs)^2)
#   sse[rr] = sum(se)
#   #lines(xshort, predict(shortfit, data.frame(x = xshort)), col = 'red')
#   }
# worst5 = order(sse)[96:100]
# preds = preds[-worst5, ]
# ssse[ss] = sum((sse))
# for(ii in 1:95){
#   lines(xshort, preds[ii,], col = tnsblack)
# }
# lines(xshort, probs, lwd = 3, col = 'blue')
#   }
# dev.off()
# plot(spans, ssse)
  
## Looks reasonable.
## I used the actual numbers of individuals of each age from Argentina, so these fits are representative of how well we might expect to approximate the true population pattern using this smoother.
## Most of the time, the soother approximates the true, underlying population relationship.
## The smoothed estimate doesn't always match the underlying population pattern perfectly, but it appears to be an unbiased estimator (I've run this several times to get a feel for the range of possible outcomes).
## The soothed estimate from the binary data usually matches the true undelrying pattern a bit better than the weighted summary smoother, so I'll use the raw, binary data to fit my model inputs
##### End exercise ##################
######################################



##  Null distributions will be similar no matter what down sampling scheme we choose
##  There are differences between observed dist of H1N1 and H3N2 from null
table(null_set_1$country, null_set_1$season)
table(null_set_2$country, null_set_2$season)
## We included the same combo of countries and seasons in both null sets. These are the countries and seasons from which we will draw data.
ccs = rep(unique(null_set_1$country), each = length(unique(null_set_1$season)))
ssns = rep(unique(null_set_1$season), times = length(unique(null_set_1$country)))
rns = paste(ccs, ssns, sep = "_") # Row names of the master matrix



## Make a master table of observed case counts
H1.master = H3.master = matrix(0, nrow = length(rns), ncol = length(18:90), dimnames = list(rns, 18:90))
for(rr in 1:length(rns)){
  # Extract data from country-season of interest as a data frame
  valid = subset(dat.002, country == ccs[rr] & season == ssns[rr])
  # Count up H1N1 cases in each age group from that country and season
  H1.master[rr, ] = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$flutype == 1)})
  # Count up H3N2 cases in each age group from that country and season
  H3.master[rr, ] = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$flutype == 2)})
}

## Use smoothed density to form null predictions
null.master1 = null.master2 = H1.master*0
for(rr in 1:length(rns)){
  # Extract data from country of interest as a data frame
  valid1 = subset(null_set_1, country == ccs[rr])
  # Calculate smoothed density using default bandwidth
  ad1 = density(x = valid1$age, from = 18, to = 90, n = 73)
  ad1 = ad1$y/sum(ad1$y)
  ad1 = matrix(ad1, nrow = length(unique(ssns)), ncol = 73, byrow = T)
  null.master1[ccs == ccs[rr], ] = ad1
  ## Repeat for 2nd null set
  valid2 = subset(null_set_2, country == ccs[rr])
  # Calculate smoothed density using default bandwidth
  ad2 = density(x = valid2$age, from = 18, to = 90, n = 73)
  ad2 = ad2$y/sum(ad2$y)
  ad2 = matrix(ad2, nrow = length(unique(ssns)), ncol = 73, byrow = T)
  null.master2[ccs == ccs[rr], ] = ad2
}
barplot(null.master1)
barplot(null.master2)

## Create av, dx and vac matrices
av.master = vac.master = dx.master = H1.master*0
for(rr in 1:length(unique(ccs))){
  # Extract data from country of interest as a data frame
  valid = subset(dat.002, country == unique(ccs)[rr])
  # loess smoother using degree 2 and span = .9. Do sensitivity analyses to check these later.
  if(unique(ccs)[rr] == 'Belgium'){
  vacmod = loess(formula = anyvac ~ age, data = valid, degree = 2, span = .9)
  }else{
  vacmod = loess(formula = anyvac ~ age, data = valid, degree = 2, span = .9)
  }
  vac = predict(vacmod, data.frame(age = 18:90))
  vac = matrix(vac, nrow = length(unique(ssns)), ncol = 73, byrow = T)
  vac.master[ccs == unique(ccs)[rr], ] = vac
  ## Check that splines are sensible
  denom = sapply(18:90, FUN = function(aa){sum(valid$age == aa)})
  num = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$anyvac == 1)})
  plot(18:90, num/denom, xlab = 'age', ylab = 'frac vaccinated', ylim = c(0,1), main = unique(ccs)[rr])
  lines(18:90, vac[1,], col = 'blue')
  legend('topleft', c('data', 'loess smoother'), col = c('black', 'blue'), lty = c(NA, 1), pch = c(1, NA))  
  
# Repeat for antivirals  
  avmod = loess(formula = anyav ~ age, data = valid, degree = 2, span = .9)
  av = predict(avmod, data.frame(age = 18:90))
  av = matrix(av, nrow = length(unique(ssns)), ncol = 73, byrow = T)
  av.master[ccs == unique(ccs)[rr], ] = av
  ## Check that splines are sensible
  denom = sapply(18:90, FUN = function(aa){sum(valid$age == aa)})
  num = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$anyav == 1)})
  plot(18:90, num/denom, xlab = 'age', ylab = 'frac av use', ylim = c(0,1), main = unique(ccs)[rr])
  lines(18:90, av[1,], col = 'blue')
  legend('topleft', c('data', 'loess smoother'), col = c('black', 'blue'), lty = c(NA, 1), pch = c(1, NA))  
  
  
  # Repeat for underlying symptoms  
  dxmod = loess(formula = anydx ~ age, data = valid, degree = 2, span = .9)
  dx = predict(dxmod, data.frame(age = 18:90))
  dx = matrix(dx, nrow = length(unique(ssns)), ncol = 73, byrow = T)
  dx.master[ccs == unique(ccs)[rr], ] = dx
  ## Check that splines are sensible
  denom = sapply(18:90, FUN = function(aa){sum(valid$age == aa)})
  num = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$anydx == 1)})
  plot(18:90, num/denom, xlab = 'age', ylab = 'frac w. underlying symp', ylim = c(0,1), main = unique(ccs)[rr])
  lines(18:90, dx[1,], col = 'blue')
  legend('topleft', c('data', 'loess smoother'), col = c('black', 'blue'), lty = c(NA, 1), pch = c(1, NA))  
}

## Correct Belgium estimates (odd estimates in lower age groups where few patients enrolled)
## Strategy: Set to biologically reasonable boundaries, and don't worry too much. 
## Age-specific null predicts very few cases here, and few data points observed, so these 
## old Belgian cases won't be too influential.
av.master[which(av.master < 0)] = 0
vac.master[which(vac.master > .75)] = .75 # Assume no more than 75% vaccine coverage in any age group
# Set about replacing NAs
replaceNA = function(row){
  last = element = length(row)
  # Move backward from the end of the row until you hit the first number
  while(is.na(row[element])){
    element = element - 1
  }
  # Repeat the last estimated number for all the older age groups
  row[element:last] = row[element]
  row # return the new row
}
# av
edit = which(apply(X = av.master, MARGIN = 1, FUN = function(xx){any(is.na(xx))}))
for(rr in edit){
  av.master[rr, ] = replaceNA(av.master[rr,])
}
# vac
edit = which(apply(X = vac.master, MARGIN = 1, FUN = function(xx){any(is.na(xx))}))
for(rr in edit){
  vac.master[rr, ] = replaceNA(vac.master[rr,])
}
# dx
edit = which(apply(X = dx.master, MARGIN = 1, FUN = function(xx){any(is.na(xx))}))
for(rr in edit){
  dx.master[rr, ] = replaceNA(dx.master[rr,])
}




## Create expected protection matrices
load('INSIGHT002_weights.RData')
# For a SH season, p(by = y - 1) = 0.5
#                  p(by = y) = 0.5
# For a NH season, p(by = y1 - 1) = 0.0625
#                  p(by = y1) = 0.875
#                  p(by = y2) = 0.0625
## Extract y1 from season rownames
y1s = as.numeric(gsub(pattern = "(\\w+_\\w\\w\\.)(\\d\\d)(\\.?\\d?\\d?)", replacement = "\\2", x = rownames(null.master2)))
# Extract hemisphere
hemisphere = gsub(pattern = "(\\w+)_(\\w\\w)(\\.\\d\\d)(\\.?\\d?\\d?)", replacement = "\\2", x = rownames(null.master2))
# Initialize
proH1.master = proH2.master = proH3.master = prog1.master = prog2.master = proN1.master = proN2.master = vac.master*0

for(rr in 1:nrow(null.master2)){
  ages = as.numeric(colnames(null.master2))
  if(hemisphere[rr] == "SH"){
    b2 = y1s[rr]+2000-ages
    b1 = b2-1
    rn = paste(y1s[rr]+2000, ccs[rr], sep = "") # Extract the weights row name for y-1 and country
    #                 #- p_H1N1 from year and country -#           #- p_H1N1 from year and country -#       
    #                 #- cols represent by1           -#           #- cols represent by2           -#
  proH1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.5 + weights[[1]][rn, as.character(b2)]*0.5
  proH2.master[rr,] = weights[[2]][rn, as.character(b1)]*0.5 + weights[[2]][rn, as.character(b2)]*0.5
  proH3.master[rr,] = weights[[3]][rn, as.character(b1)]*0.5 + weights[[3]][rn, as.character(b2)]*0.5
  prog1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.5 + weights[[1]][rn, as.character(b2)]*0.5 + weights[[2]][rn, as.character(b1)]*0.5 + weights[[2]][rn, as.character(b2)]*0.5
  prog2.master[rr,] = weights[[3]][rn, as.character(b1)]*0.5 + weights[[3]][rn, as.character(b2)]*0.5
  proN1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.5 + weights[[1]][rn, as.character(b2)]*0.5
  proN2.master[rr,] = weights[[3]][rn, as.character(b1)]*0.5 + weights[[3]][rn, as.character(b2)]*0.5 + weights[[2]][rn, as.character(b1)]*0.5 + weights[[2]][rn, as.character(b2)]*0.5
  }else if(hemisphere[rr]=='NH'){
    b2 = y1s[rr]+2000-ages
    b1 = b2-1
    b3 = b2+1
    rn = paste(y1s[rr]+2000, ccs[rr], sep = "") # Extract the row name for y-1 and country
    #                 #- p_H1N1 from year and country -#           #- p_H1N1 from year and country -#       
    #                 #- cols represent by1           -#           #- cols represent by2           -#
    proH1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.0625 + weights[[1]][rn, as.character(b2)]*0.875 + weights[[1]][rn, as.character(b3)]*0.0625
    proH2.master[rr,] = weights[[2]][rn, as.character(b1)]*0.0625 + weights[[2]][rn, as.character(b2)]*0.875 + weights[[2]][rn, as.character(b3)]*0.0625
    proH3.master[rr,] = weights[[3]][rn, as.character(b1)]*0.0625 + weights[[3]][rn, as.character(b2)]*0.875 + weights[[3]][rn, as.character(b3)]*0.0625
    
    prog1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.0625 + weights[[1]][rn, as.character(b2)]*0.875 + weights[[1]][rn, as.character(b3)]*0.0625 + weights[[2]][rn, as.character(b1)]*0.0625 + weights[[2]][rn, as.character(b2)]*0.875 + weights[[2]][rn, as.character(b3)]*0.0625
    prog2.master[rr,] = weights[[3]][rn, as.character(b1)]*0.0625 + weights[[3]][rn, as.character(b2)]*0.875 + weights[[3]][rn, as.character(b3)]*0.0625
    proN1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.0625 + weights[[1]][rn, as.character(b2)]*0.875 + weights[[1]][rn, as.character(b3)]*0.0625
    proN2.master[rr,] = weights[[3]][rn, as.character(b1)]*0.0625 + weights[[3]][rn, as.character(b2)]*0.875 + weights[[3]][rn, as.character(b3)]*0.0625 + weights[[2]][rn, as.character(b1)]*0.0625 + weights[[2]][rn, as.character(b2)]*0.875 + weights[[2]][rn, as.character(b3)]*0.0625
  }
}


par(mfrow = c(3,1))
barplot(null.master2, main = "All tested")
barplot(H1.master, main = "confirmed H1")
barplot(H3.master, main = "confirmed H3")
barplot(av.master, main = 'antiviral use')
barplot(dx.master, main = 'underlying symptoms present')
barplot(vac.master, main = 'vaccinated')
par(mfrow = c(2,1))
barplot(proH1.master, main = "H1 protection")
barplot(proH3.master, main = "H3 protection")
barplot(prog1.master, main = "g1 protection")
barplot(prog2.master, main = "g2 protection")
barplot(proN1.master, main = "N1 protection")
barplot(proN2.master, main = "N2 protection")

