setwd('~/Dropbox/R/2017_INSIGHT/')
rm(list = ls())
# Load data
dat.003 = read.csv('data003.csv', colClasses = c('character', 'integer', 'integer', 'integer', 'character', 
                                                 'integer', 'character', 'character', 'character', 'character', 
                                                 'character', 'integer', 'character', 'integer', 'integer', 
                                                 'integer', 'integer'))


# Enter actual country names instead of country codes
Country = character(nrow(dat.003))
ccs = unique(dat.003$COUNTRY_CODE); ccs # Existing codes
cnames = c('Denmark', 'Spain', 'Germany', 'USA', 'Belgium', 'Norway', 'Poland', 'Austria', 'UK', 'Australia', 'China', 'Thailand', 'Argentina', 'Chile', 'Greece', 'Peru', 'Singapore') # Names with which to replace 
rbind(ccs, cnames) # Check
for(ii in 1:length(ccs)){
  Country[which(dat.003$COUNTRY_CODE == ccs[ii])] = cnames[ii]
}
# # Check
#sm = sample(x = 1:nrow(dat.003), size = 150, replace = FALSE)
# rbind(dat.003$COUNTRY_CODE[sm], Country[sm])
dat.003$country = Country
rm(cnames, ccs, Country)

## reformat year:
dat.003$year = dat.003$year+2000 # Change from 2 didigit to 4 digit

## Estimate two possible birth years and probabilities of birth in that year, given month of case observation. Add this info to data frame.
dat.003$by1 = dat.003$year - dat.003$age - 1 # First possible birth year
dat.003$by2 = dat.003$year - dat.003$age  #Second possible birth year

# Exclude cases born before 1918
dat.003 = dat.003[-which(dat.003$by1 < 1918), ]
dat.003 = dat.003[-which(is.na(dat.003$by1)), ]

# Input probabilities of birth in candidate year 1 or 2 given month of case observation:
p.by1.given.month = seq(11.5, .5, -1)/12
p.by2.given.month = seq(.5, 11.5, 1)/12

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
# Make sure all month entries follow the specified format
if(any(!(dat.003$month %in% c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov','Dec')))){ error('Ivalid month passed to month_to_num')}
# Convert
dat.003$num_month = month_to_num(dat.003$month)

# Calculate data frame values
dat.003$p.by1 = p.by1.given.month[dat.003$num_month]
dat.003$p.by2 = p.by2.given.month[dat.003$num_month]



# Calculate imprinting probs based on birth year
table(dat.003$COUNTRY_CODE, dat.003$year) # View countires and years in which we need to do reconstrutions
# source('../Reconstructions/Infection_age_structure_2017.R')
# weights = get.type.weights.AB(years.out = 2009:2017, Countries.out =  c('Denmark', 'Spain', 'Germany', 'USA', 'Belgium', 'Norway', 'Poland', 'Austria', 'UK', 'Australia', 'China', 'Thailand', 'Argentina', 'Chile', 'Greece', 'Peru', 'Singapore'), region.in = c('Euro', 'Euro', 'Euro', 'default', 'Euro', 'Euro', 'Euro', 'Euro', 'Euro', 'default', 'Asia', 'Asia', 'default', 'default', 'Euro', 'default', 'Asia'))
# save(weights, file = 'INSIGHT003_weights.RData')
load('INSIGHT003_weights.RData')


# Assign individuals a probability of imprinting protection against g1 or g2
# Calculate probabilities as a weighted average of p(x|by1) and p(x|by2)
rnms = paste(dat.003$year, dat.003$country, sep = '')
dat.003$p.g1.protection = mapply(function(rn, cn1, cn2, p1, p2){ (weights[[1]][rn, cn1]+weights[[2]][rn,cn1])*p1 + (weights[[1]][rn, cn2]+weights[[2]][rn,cn2])*p2 }, rn = rnms, cn1 = as.character(dat.003$by1), cn2 = as.character(dat.003$by2), p1 = dat.003$p.by1, p2 = dat.003$p.by2 ) 

dat.003$p.g2.protection = mapply(function(rn, cn1, cn2, p1, p2){ (weights[[3]][rn, cn1])*p1 + (weights[[3]][rn, cn2])*p2 }, rn = rnms, cn1 = as.character(dat.003$by1), cn2 = as.character(dat.003$by2), p1 = dat.003$p.by1, p2 = dat.003$p.by2 ) 

# Assign individuals a probability of imprinting protection against H.subtype
dat.003$p.H1.protection = mapply(function(rn, cn1, cn2, p1, p2){ (weights[[1]][rn, cn1])*p1 + (weights[[1]][rn, cn2])*p2 }, rn = rnms, cn1 = as.character(dat.003$by1), cn2 = as.character(dat.003$by2), p1 = dat.003$p.by1, p2 = dat.003$p.by2 ) 
dat.003$p.H3.protection = mapply(function(rn, cn1, cn2, p1, p2){ (weights[[3]][rn, cn1])*p1 + (weights[[3]][rn, cn2])*p2 }, rn = rnms, cn1 = as.character(dat.003$by1), cn2 = as.character(dat.003$by2), p1 = dat.003$p.by1, p2 = dat.003$p.by2 ) 

# Assign individuals a probability of imprinting protection against N.sub
dat.003$p.N1.protection = mapply(function(rn, cn1, cn2, p1, p2){ (weights[[1]][rn, cn1])*p1 + (weights[[1]][rn, cn2])*p2 }, rn = rnms, cn1 = as.character(dat.003$by1), cn2 = as.character(dat.003$by2), p1 = dat.003$p.by1, p2 = dat.003$p.by2 ) 
dat.003$p.N2.protection = mapply(function(rn, cn1, cn2, p1, p2){ (weights[[3]][rn, cn1]+weights[[2]][rn, cn1])*p1 + (weights[[3]][rn, cn2]+weights[[2]][rn, cn2])*p2 }, rn = rnms, cn1 = as.character(dat.003$by1), cn2 = as.character(dat.003$by2), p1 = dat.003$p.by1, p2 = dat.003$p.by2 ) 

# Create a categorical variable that tracks high, med and low probs of protection
dat.003$g1.category = dat.003$g2.category = 'med'
dat.003$g1.category[which(dat.003$p.g1.protection <.1)] = 'low'
dat.003$g1.category[which(dat.003$p.g1.protection >.9)] = 'high'
dat.003$g2.category[which(dat.003$p.g2.protection <.1)] = 'low'
dat.003$g2.category[which(dat.003$p.g2.protection >.9)] = 'high'






# Define seasons based on October to October years
sort(unique(dat.003$year))
NHfirst.part = c('Oct', 'Nov', 'Dec')
NHsecond.part = c('Jan', 'Feb', 'Mar')
SH = c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
ssn = vector('character', nrow(dat.003)) # Create a new vector to record season codes
season.codes.NH = paste('NH', c('09', 10:16), 10:17, sep = '.')
season.codes.SH = paste('SH', 10:17, sep = '.')
# Assign each observation to a NH or SH season (Oct-March = NH, April - Sept = SH)
yy = 2009
for(ii in 1:8){
  NH.indices = c(which(dat.003$year == yy & dat.003$month %in% NHfirst.part), which(dat.003$year == yy+1 &  dat.003$month %in% NHsecond.part))
  ssn[NH.indices] = season.codes.NH[ii]
  SH.indices = which(dat.003$year == yy + 1 & dat.003$month %in% SH)
  ssn[SH.indices] = season.codes.SH[ii]
  yy = yy + 1
}
# Add the season assignements to the data frame
dat.003$season = (ssn)
rm(ssn)


# Create a categorical variable that tracks high, med and low probs of protection
dat.003$g1.category = dat.003$g2.category = 'med'
dat.003$g1.category[which(dat.003$p.g1.protection <.1)] = 'low'
dat.003$g1.category[which(dat.003$p.g1.protection >.9)] = 'high'
dat.003$g2.category[which(dat.003$p.g2.protection <.1)] = 'low'
dat.003$g2.category[which(dat.003$p.g2.protection >.9)] = 'high'








## Meanwhile, also calculate the probability that each case would have been challenged by H1N1 vs. H3N2, given the country and year in which their case was observed
history = read.csv(file = 'Historical_circulation_data.csv', stringsAsFactors = FALSE)
# STore historical circulation data in this matrix
hist_H1fracs = matrix(NA, nrow = length(season.codes.NH)+length(season.codes.SH), ncol = length(unique(dat.003$country)), dimnames = list(NULL, unique(dat.003$country)))
rownns = vector('numeric', nrow(hist_H1fracs))
# Populate the matrix
yy = 2010
yy2 = 10
rr = 1
for(ii in 1:8){
  # Extract overall NH and SH data for seasons with insufficient country-specific data
  valid.overall.NH = subset(history, (Year == yy-1 & Week >=40)|(Year == yy & Week <14))
  valid.overall.SH = subset(history, (Year == yy & Week >= 14 & Week < 40))
  for(cc in 1:length(unique(dat.003$country))){
    #Extract surveillance data from weeks belonging to the NH season ending in year yy
    valid.NH = subset(history, (Year == yy-1 & Week >=40 & Country == colnames(hist_H1fracs)[cc])|(Year == yy & Week <14 & Country == colnames(hist_H1fracs)[cc]))
    # If more than 30 observations, calculate and store the H1 frac
    if(dim(valid.NH)[1]>0){
    if(sum(valid.NH[,4:6], na.rm = TRUE)>=30){
      hist_H1fracs[rr, cc] = sum(valid.NH[,4:5], na.rm = TRUE)/sum(valid.NH[,4:6], na.rm = TRUE)
    }else{ # Else, fill in the overall measurement
      hist_H1fracs[rr, cc] = sum(valid.overall.NH[,4:5], na.rm = TRUE)/sum(valid.overall.NH[,4:6], na.rm = TRUE)
    }
    }else{ # Else, fill in the overall measurement
      hist_H1fracs[rr, cc] = sum(valid.overall.NH[,4:5], na.rm = TRUE)/sum(valid.overall.NH[,4:6], na.rm = TRUE)
    }
    valid.SH = subset(history, (Year == yy & Week >= 14 & Week < 40 & Country == colnames(hist_H1fracs)[cc]))
    # Calculate and store the H1 frac
    if(dim(valid.SH)[1]>0){
    if(sum(valid.SH[,4:6], na.rm = TRUE)>=30){
      hist_H1fracs[rr+1, cc] = sum(valid.SH[,4:5], na.rm = TRUE)/sum(valid.SH[,4:6], na.rm = TRUE)
    }else{
      hist_H1fracs[rr+1, cc] = sum(valid.overall.SH[,4:5], na.rm = TRUE)/sum(valid.overall.SH[,4:6], na.rm = TRUE)
    }
    }else{
      hist_H1fracs[rr+1, cc] = sum(valid.overall.SH[,4:5], na.rm = TRUE)/sum(valid.overall.SH[,4:6], na.rm = TRUE)
    } # End if statement
  } # End country loop
  rownns[rr] = paste('NH', yy2-1, yy2, sep = '.')
  rownns[rr+1] = paste('SH', yy2, sep = '.')
  rr = rr + 2
  yy = yy + 1
  yy2 = yy2+1
} # End year loop
rownns[1] = "NH.09.10"
rownames(hist_H1fracs)=rownns


## EXCLUDE CASES OF:
# FLUB
# NOT SUBTYPED
# COINFECTION
# NA
unique(dat.003$flutype)
exclude = which(dat.003$flutype %in% c(3,6,NA))
dat.003 = dat.003[-exclude, ]



## Assign challenges to each case
challenge = vector('character', nrow(dat.003))
challenge[which(dat.003$flutype==1)] = 'H1N1' # If confirmed H1N1, definitely H1N1 challenge
challenge[which(dat.003$flutype==2)] = "H3N2" # If confirmed H3N2, definitely H3N2 challenge
# Fill in confirmed PCR negative probabilistically:
fillin = which(dat.003$flutype==4)
for(ii in fillin){
  pH1 = hist_H1fracs[dat.003$season[ii], dat.003$country[ii]] # Prob H1 = fraction of all circulating A viruses of H1N1 subtype
  challenge[ii] = sample(x = c("H1N1", 'H3N2'), size = 1, prob = c(pH1, 1-pH1)) # Randomly draw a challenge using sample prob
}
dat.003$challenge = challenge



## Add column for H1N1 (Confirmed H1N1 only)
dat.003$H1N1 = 0; dat.003$H1N1[which(dat.003$flutype == 1)] = 1; dat.003$H1N1 = as.numeric(dat.003$H1N1)
## Add column for H3N2 (Confirmed H3N2 only)
dat.003$H3N2 = 0; dat.003$H3N2[which(dat.003$flutype == 2)] = 1; dat.003$H3N2 = as.numeric(dat.003$H3N2)
## Add column for Inefted (Confirmed H1N1 or H3N2)
dat.003$infected = 0; dat.003$infected[which(dat.003$flutype %in% c(1,2))] = 1; dat.003$infected = as.numeric(dat.003$infected)
## Add column for protected_grp (p(protection @ group level) given imputed challenge)
dat.003$protected_grp = dat.003$p.g1.protection # First, fill in g1 protection probs
overwrite = which(dat.003$challenge == 'H3N2') # STore indices of H3N2 challenges 
dat.003$protected_grp[overwrite] = dat.003$p.g2.protection[overwrite] # Overwrite with g2 probs
## Add column for protected_sub (p(protection @ subtype level) given imputed challenge)
dat.003$protected_sub = 0
dat.003$protected_sub[which(dat.003$challenge == 'H1N1')] = dat.003$p.H1.protection[which(dat.003$challenge == 'H1N1')]
dat.003$protected_sub[which(dat.003$challenge == 'H3N2')] = dat.003$p.H3.protection[which(dat.003$challenge == 'H3N2')]
## Add column for protected_N (p(protection @ NA grp level) given imputed challenge)
dat.003$protected_N = 0
dat.003$protected_N[which(dat.003$challenge == 'H1N1')] = dat.003$p.N1.protection[which(dat.003$challenge == 'H1N1')]
dat.003$protected_N[which(dat.003$challenge == 'H3N2')] = dat.003$p.N2.protection[which(dat.003$challenge == 'H3N2')]


write.table(x = dat.003, file = "Clean_dat003.txt")

