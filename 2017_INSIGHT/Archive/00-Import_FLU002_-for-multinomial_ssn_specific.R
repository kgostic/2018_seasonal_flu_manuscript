## This file loads INSIGHT FLU002 data and inputs for multinomial analysis.
##   Only analyze seasons in which >= 50 cases of each subtype were confirmed
##   Fit age-specific risk curves specific to each season


###############################
# 1. Load raw data
###############################
dat.002 = read.csv('raw-data/data002.csv', colClasses = c('character', 'integer', 'integer', 'integer', 'character', 'integer', 'character', 'character', 'character', 'character', 'character', 'character', 'character', 'character'))
head(dat.002)



###############################
# 2. Clean raw data
###############################
Country = character(nrow(dat.002))
ccs = unique(dat.002$COUNTRY_CODE); ccs # Existing codes
cnames = c('Denmark', 'Spain', 'Germany', 'Estonia', 'USA', 'Belgium', 'Portugal', 'Poland', 'Austria', 'UK', 'Australia', 'Thailand', 'Argentina', 'Chile', 'Greece', 'Peru', 'Japan') # Names with which to replace 
rbind(ccs, cnames) # Check
# Add true country names to data frame
for(ii in 1:length(ccs)){
  Country[which(dat.002$COUNTRY_CODE == ccs[ii])] = cnames[ii]
}
dat.002$country = Country
rm(cnames, ccs, Country) # Remove code names
# Exclude cases over age 90
# Do this because individuals observed in 2009 of age 90 could have been born in 1918
# Older individuals could have been born before 1918
paste('Excluding', sum(dat.002$age > 90, na.rm = TRUE), 'cases over age 90 from INSIGHT 002 data.')
dat.002 = dat.002[-which(dat.002$age > 90), ]

# Make sure all month entries follow the specified format
if(any(!(dat.002$month %in% c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov','Dec')))){ error('Ivalid month passed to month_to_num')}
# Input numeric months
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

# Define seasons based on October to October years
sort(unique(dat.002$year))
NHfirst.part = c('Oct', 'Nov', 'Dec')
NHsecond.part = c('Jan', 'Feb', 'Mar')
SH = c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')
ssn = vector('character', nrow(dat.002)) # Initialize an empty vector to record season codes
season.codes.NH = paste('NH', c('09', 10:16), 10:17, sep = '.') # Genearte character vectors of NH and SH seasons in data
season.codes.SH = paste('SH', 10:17, sep = '.')
# Assign each observation in the dat.002 data frame to a NH or SH season (Oct-March = NH, April - Sept = SH)
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
rm(ssn) # Tidy up

## Pull out variables used in model
dat.002 = subset(dat.002, select = c('age', 'anyvac', 'anydx', 'anyav', 'season', 'country', 'season', 'flutype'))
nrow(dat.002) -> original

## Drop rows with missing data:
paste('Dropping', sum((apply(dat.002, MARGIN = 1,  function(xx) (any(is.na(xx)))))), 'rows with NA')
# Exclude rows with any NAs
dat.002 = na.omit(dat.002)
# Exclude those with unknown vaccination status.
invalid = which(dat.002$anyvac == 2)
dat.002 = dat.002[-invalid, ]
paste('Dropping', length(invalid), 'cases wtih unknown vaccination status.')
# Exclude two observations with very low age
invalid = which(dat.002$age == 17)
dat.002 = dat.002[-invalid,]
paste('Dropping', length(invalid), 'cases under age 18.')
## Exclude cases with unknown flutype or coinfection
dat.002 = subset(dat.002, flutype %in% c(1,2,3,4))
paste('Excluded', original-nrow(dat.002), 'total cases with missing info from dat.002')




#### EXCLUDE ALL OBSERVATIONS, EXCEPT IN SEASONS IN WHICH BOTH SUBTYPES CIRCULATED
tb = table(dat.002$season, dat.002$flutype)
keep = which(tb[,1]>=50 & tb[,2]>=50) # Indices of rows in tb that meet the critera
keepssns = rownames(tb[keep,]) # Extract names of valid seasons
dat.002 = subset(dat.002, season %in% keepssns) # Drop cases from all other seasons

## Check data
# check case counts by country and subtype
table(dat.002$country, dat.002$flutype) -> tab1
tab1
tab1[,1:2]
rowSums(tab1[,1:2])

############# END DATA CLEANING #######################
############# dat.002 contains clean data #############







##############################################
##  4. Reformat data for input into likelihood
##     Generate corresponding inputs for each tested factor:
##      - Vaccination in the past year
##      - Antiviral treatment
##      - Presence of underlying conditions
##      - Age (fit as a step function)
##      - Prob of imprinting at the HA group level
##      - " ......................" HA subtype level
##      - " ......................" NA subtype level
##############################################
## rns lists all unique combinations of season and country
## rns will form row names in a master matrix where we store data and relevant model inputs for each country-season of case observation
ccs = rep(unique(dat.002$country), each = length(unique(dat.002$season)))
ssns = rep(unique(dat.002$season), times = length(unique(dat.002$country)))
rns = paste(ccs, ssns, sep = "_") # Row names of the master matrix
####################################################################
## 4a. Make a master table of observed case counts in each country, season and single year of age. Country-years go on rows, age goes on columns.
##  H1.master - confirmed H1N1 case counts
##  H3.master - confirmed H3N2 case counts
##  tested.master - all subjects tested (counts)
##  av.master - fraction of tested cases in each country-season, and single year of age that was treated with antivirals
##  dx.master - fraction of tested cases in each country-season and single year of age that had underlying conditions
##  vac.master - fraction of tested cases in each country-season and single year of age that were vaccinated
H1.master = H3.master = tested.master = av.master = vac.master =  dx.master = matrix(0, nrow = length(rns), ncol = length(18:90), dimnames = list(rns, 18:90))
for(rr in 1:length(rns)){
  # Extract data from country-season of interest as a data frame
  valid = subset(dat.002, country == ccs[rr] & season == ssns[rr])
  # Count the number of subjects from the focal country-season in each single year of age from 18-90
  tested = sapply(18:90, FUN = function(aa){sum(valid$age == aa)})
  tested.master[rr, ] = tested ## Store the country-year specific counts into the master matrix
  
  ## Fill in av.master
  num = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$anyav == 1)}) # Count the number of subjects in each age bin treated with antivirals
  av.master[rr, ] = (num/ifelse(tested > 0, tested, 1)) # Divide by the number tested to get the proprtion treated with avs in each age bin. 
## Note - if no patients were enrolled in a particular age bin, normalizing leads to a divide by zero error.
## To avoid this, substitute a 1 where tested == 0 to avoid 0 denominator.
## This substitituion will only be made when no patients were enrolled in a particular age group, so the numerator (number treated with antivirals) will also be 0. The final entry in av.master will be 0/1 = 0
  
  ## Repeat for proportion vaccinated:
  num = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$anyvac == 1)}) ## Fill in master matrix
  vac.master[rr, ] = (num/ifelse(tested > 0, tested, 1))
  
  ## Repeat for proportion with underlying conditions:
  num = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$anydx == 1)})
  dx.master[rr, ] = (num/ifelse(tested > 0, tested, 1))
  
  # Count up H1N1 cases in each age group from that country and season
  H1.master[rr, ] = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$flutype == 1)})
  # Count up H3N2 cases in each age group from that country and season
  H3.master[rr, ] = sapply(18:90, FUN = function(aa){sum(valid$age == aa & valid$flutype == 2)})
}


##############################################
##  4b. Import imprinting reconstructions
## (i.e. probs of imprinting to a particular subtype, given birth year)
##############################################
table(dat.002$COUNTRY_CODE, dat.002$year) # View countires and years in which we need to do reconstrutions
## See ../Reconstructions/02-create-INSIGHT-weights.R for code to do reconstructions
load('../../Reconstructions/INSIGHT_weights2018-12-04.RData')
## this loads a variable named weights
## weights is a list with elements:
## weights$weights.master.1 -> probabilities of imprinting to H1N1 in childhood
## weights$weights.master.2 -> probabilities of imprinting to H2N2 in childhood
## weights$weights.master.3 -> probabilities of imprinting to H3N2 in childhood
## weights$weights.master.naiive -> probabilities of not yet having imprinted (nonzero probs in children only)



####################################################################
## 4c. create corresponding tables of the probability of imprinting protection 
##     for individuals of age a, in a given country-season:
##       - First, determine possible birth years based on age at the time of study enrollment.
##       - Then, calculate probabilities of protection. p(protection|birth year) are found in the weights list loaded above.
################
##  Rationale for multiple possible birth years:
# SH season includes months April-Sept
##  Given age observed April-Sept, there's a 50% chance that your birthday was in the previous calendar year, and a 50% chance that your birthday was in the current calendar year.
##     --> A newborn of age 0 in April could have been born in May-Dec. of the previous year (8/12 chance of being born in year y-1), or could have been born in Jan-April of the curren tyear (4/12 chance of being born in year y)
##     --> Following the same logic, and averaging across probabilities of birth in year y, or y-1, given observation in April-Sept, the overall probabilities of birth in years y and y-1 come out to:
##                 p(by = y - 1) = 0.5
#                  p(by = y) = 0.5
# For a NH season, your age is observed at the end of the calendar, year, so there's a lower chance that your birthday was in the previous calendar year, and a higher chance that your birthday is in the current calendar year.
##  Follow the same logic as above
##  See Methods for a better description
##                 p(by = y1 - 1) = 0.0625
#                  p(by = y1) = 0.875
#                  p(by = y2) = 0.0625
## Extract y1 from season rownames.  (y1 = year of case observation, or first year relevant to a NH season. E.g. y1 = 2000 for the 2000 SH season, or for the 2000-2001 NH season)
y1s = as.numeric(gsub(pattern = "(\\w+_\\w\\w\\.)(\\d\\d)(\\.?\\d?\\d?)", replacement = "\\2", x = rownames(H1.master)))
# Extract hemisphere
hemisphere = gsub(pattern = "(\\w+)_(\\w\\w)(\\.\\d\\d)(\\.?\\d?\\d?)", replacement = "\\2", x = rownames(H1.master))
# Initialize matrices of the same dimensions as the data master matrices to store imprinting protection probs
proH1.master = proH2.master = proH3.master = prog1.master = prog2.master = proN1.master = proN2.master = H1.master*0

## For each country-season of case observation (i.e. each row of the master matrix),
##  calculate age-specific probs of imprinting protection
for(rr in 1:nrow(H1.master)){
  ages = as.numeric(colnames(H1.master)) # Get relevant single years of age
  if(hemisphere[rr] == "SH"){
    b2 = y1s[rr]+2000-ages # Later possible birth year given by current year-age
    b1 = b2-1 ## Also possible to be born one year before y1-age (e.g. if your birthday is late in the year)
    rn = paste(y1s[rr]+2000, ccs[rr], sep = "") # Extract the weights row name for y-1 and country
    #                 #- p_H1N1 from year and country -#           #- p_H1N1 from year and country -#       
    #                 #- cols represent by1           -#           #- cols represent by2           -#
  proH1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.5 + weights[[1]][rn, as.character(b2)]*0.5 ## HA subtype-specific protection against H1
  proH2.master[rr,] = weights[[2]][rn, as.character(b1)]*0.5 + weights[[2]][rn, as.character(b2)]*0.5 ## H2 subtype-specific protection
  proH3.master[rr,] = weights[[3]][rn, as.character(b1)]*0.5 + weights[[3]][rn, as.character(b2)]*0.5 ## HA subtype-specific protection against H3
  prog1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.5 + weights[[1]][rn, as.character(b2)]*0.5 + weights[[2]][rn, as.character(b1)]*0.5 + weights[[2]][rn, as.character(b2)]*0.5 ## HA group-specific protection against H1
  prog2.master[rr,] = weights[[3]][rn, as.character(b1)]*0.5 + weights[[3]][rn, as.character(b2)]*0.5 ## HA group-specific protection against H3
  proN1.master[rr,] = weights[[1]][rn, as.character(b1)]*0.5 + weights[[1]][rn, as.character(b2)]*0.5 ## NA subtype-specific protection against N1
  proN2.master[rr,] = weights[[3]][rn, as.character(b1)]*0.5 + weights[[3]][rn, as.character(b2)]*0.5 + weights[[2]][rn, as.character(b1)]*0.5 + weights[[2]][rn, as.character(b2)]*0.5 ## NA subtype-specific protection against N2
  
  ## For NH, implement different weighting of possible birth years
  }else if(hemisphere[rr]=='NH'){
    b2 = y1s[rr]+2000-ages # Second possible birth year given by current year-age
    b1 = b2-1 # First possible birth year is in the previous year
    b3 = b2+1 # Third possible birth year is in the following year (e.g. if it's the 2000-2001 NH season, a child of age 0 could have been born in December 1999 (if the case was observed from Oct-Nov), born December 2000 (if the case was observed after December), or born January 2001). b1, b2 and b3 represent 3 possible birth years for each age observed.
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







####################################################################
## 5. Clean model inputs
####################################################################
totcases = rbind(rowSums(H1.master), rowSums(H3.master))
## Exlude country-seasons in which 0 cases were confirmed (all participating countries did not necessarily report cases to the INSIGHT databse in every year of the study)
exclude = which(colSums(totcases)== 0)
H1.master = H1.master[-exclude, ]
H3.master = H3.master[-exclude, ]
av.master = av.master[-exclude, ]
dx.master = dx.master[-exclude, ]
vac.master = vac.master[-exclude, ]
proH1.master = proH1.master[-exclude, ]
proH2.master = proH2.master[-exclude, ]
proH3.master = proH3.master[-exclude, ]
proN1.master = proN1.master[-exclude, ]
proN2.master = proN2.master[-exclude, ]
prog1.master = prog1.master[-exclude, ]
prog2.master = prog2.master[-exclude, ]
tested.master = tested.master[-exclude, ]
ccs = ccs[-exclude]
ssns = ssns[-exclude]
rns = rns[-exclude]



####################################################################
## 6. Create age-specific indicator matrices, which we will use to implement
##     age-specific step functions
####################################################################
#### Age groups of interest
## 18-24
## 25-31
## 32-38
## 39-45
## 46-52
## 53-59
## 60-66
## 67-73
## 74-80
## 81-90
a18.24 = a25.31 = a32.38 = a39.45 = a46.52 = a53.59 = a60.66 = a67.73 = a74.80 = a81.90 = H1.master*0 # Initialize
a18.24[,as.character(18:24)] = 1
a25.31[,as.character(25:31)] = 1
a32.38[,as.character(32:38)] = 1
a39.45[,as.character(39:45)] = 1
a46.52[,as.character(46:52)] = 1
a53.59[,as.character(53:59)] = 1
a60.66[,as.character(60:66)] = 1
a67.73[,as.character(67:73)] = 1
a74.80[,as.character(74:80)] = 1
a81.90[,as.character(81:90)] = 1
# Check that one and only one entry in each indicator matrix is set to 1
any(a18.24+a25.31+a32.38+a39.45+a46.52+a53.59+a60.66+a67.73+a74.80+a81.90 != 1)

# barplot(a18.24) # 7
# barplot(a25.31) # 
# barplot(a32.38)
# barplot(a39.45)
# barplot(a46.52)
# barplot(a53.59)
# barplot(a60.66)
# barplot(a67.73)
# barplot(a74.80)
# barplot(a81.90)
# barplot(a18.24+a25.31+a32.38+a39.45+a46.52+a53.59+a60.66+a67.73+a74.80+a81.90)

## Create season-specific indicator matrices, which each take value 1 in the season of interest and 0 otherwise
NH.12.13.master = SH.13.master = NH.13.14.master = NH.15.16.master = SH.16.master = H1.master*0
NH.12.13.master[grep('NH.12.13', rownames(H1.master)),] = 1
SH.13.master[grep('SH.13', rownames(H1.master)),] = 1
NH.13.14.master[grep('NH.13.14', rownames(H1.master)),] = 1
NH.15.16.master[grep('NH.15.16', rownames(H1.master)),] = 1
SH.16.master[grep('SH.16', rownames(H1.master)),] = 1
## Check that one and only one entry in each matrix takes value 1.
any(NH.12.13.master+NH.13.14.master+NH.15.16.master+SH.13.master+SH.16.master != 1)

## Check model inputs
# par(mfrow = c(3,1))
# barplot(H1.master, main = "confirmed H1")
# barplot(H3.master, main = "confirmed H3")
# 
# par(mfrow = c(2,1))
# barplot(proH1.master, main = "H1 protection")
# barplot(proH3.master, main = "H3 protection")
# barplot(prog1.master, main = "g1 protection")
# barplot(prog2.master, main = "g2 protection")
# barplot(proN1.master, main = "N1 protection")
# barplot(proN2.master, main = "N2 protection")




## Get mean age of H1N1 infection for each season
H1.age.mean = sapply(X = unique(dat.002$cs), FUN = function(ccss) mean(  subset(x = dat.002, cs == ccss & flutype == 1)$age, na.rm = TRUE  ))
H3.age.mean = sapply(X = unique(dat.002$cs), FUN = function(ccss) mean(  subset(x = dat.002, cs == ccss & flutype == 2)$age, na.rm = TRUE  ))
yr = as.numeric(sub(pattern = '.+(\\d\\d)', replacement = '\\1', x = unique(dat.002$cs)))+2000
plot(yr, H1.age.mean)

