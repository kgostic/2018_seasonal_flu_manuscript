#### Visualize two kinds of sequence data and test GAMMS
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')

## Load Arizona data
H1_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H1.csv')[,-1]
H3_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H3.csv')[,-1]
H1_AZ$year = H1_AZ$year + 1
H3_AZ$year = H3_AZ$year + 1


## Plot H1 and H3, color by year
base.cols = c('magenta', 'olivedrab1', 'red', 'yellow', 'blue', 'orchid', 'turquoise', 'goldenrod', 'purple', 'chocolate1', 'limegreen', 'dodgerblue')
cols = 'start'
for(ii in 1:12){
  ramp <- colorRamp(c(base.cols[ii], "white"))
  cols = c(cols, rgb( ramp(seq(0, .5, length = 2)), max = 255))
}
cols = cols[-1]

#reord = rev(c(1, 5, 9, 13, 17, 21, c(1, 5, 9, 13, 17, 21)+1, c(1, 5, 9, 13, 17, 21)+2, c(1, 5, 9, 13, 17, 21)+3))
reord = rev(c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24))
cols = cols[reord]

yrs = as.character(1990:2014)


####### Birth year
par(mfrow = c(2,2), mar = c(3, 3, 1, 1))
## H1 AZ
yys = seq(1920, 2000, by = 10)

bys = min(H1_AZ$BirthYear):2010; yrs = unique(H1_AZ$Year)
dt = matrix(NA, length(yrs), length(bys), dimnames = list(yrs, bys))
for(yy in 1:length(yrs)){
for(ii in 1:length(bys)){
  dt[yy, ii] = sum(H1_AZ$BirthYear == bys[ii] & H1_AZ$Year == yrs[yy])
}}
#dt = table(H1_AZ$Year, H1_AZ$BirthYear); 
use.cols = which(as.character(1990:2014) %in% rownames(dt))
xx = barplot(dt, col = cols[use.cols], border = NA, main = 'H1 Arizona', space = 0, xaxt = 'n')
use.ax = which(colnames(dt) %in% as.character(yys))
nms = colnames(dt)[use.ax]
axis(1, at = xx[use.ax], labels = nms)


## H1 NCBI
#dt = table(H1_NCBI$Year, H1_NCBI$BirthYear)
bys = min(H1_NCBI$BirthYear):2010; yrs = sort(unique(H1_NCBI$Year))
dt = matrix(NA, length(yrs), length(bys), dimnames = list(yrs, bys))
for(yy in 1:length(yrs)){
  for(ii in 1:length(bys)){
    dt[yy, ii] = sum(H1_NCBI$BirthYear == bys[ii] & H1_NCBI$Year == yrs[yy])
  }}
use.cols = which(as.character(1990:2014) %in% rownames(dt))
xx = barplot(dt, col = cols[use.cols], border = NA, main = 'H1 NCBI', space = 0, xlab = '', xaxt = 'n')
use.ax = which(colnames(dt) %in% as.character(yys))
nms = colnames(dt)[use.ax]
axis(1, at = xx[use.ax], labels = nms)



## H3 AZ
#dt = table(H3_AZ$Year, H3_AZ$BirthYear); use.cols = which(as.character(1990:2014) %in% rownames(dt))
bys = min(H3_AZ$BirthYear):2010; yrs = sort(unique(H3_AZ$Year))
dt = matrix(NA, length(yrs), length(bys), dimnames = list(yrs, bys))
for(yy in 1:length(yrs)){
  for(ii in 1:length(bys)){
    dt[yy, ii] = sum(H3_AZ$BirthYear == bys[ii] & H3_AZ$Year == yrs[yy])
  }}
use.cols = which(as.character(1990:2014) %in% rownames(dt))
xx = barplot(dt, col = cols[use.cols], border = NA, main = 'H3 Arizona', space = 0, xaxt = 'n')
use.ax = which(colnames(dt) %in% as.character(yys))
nms = colnames(dt)[use.ax]
axis(1, at = xx[use.ax], labels = nms)

## H3 NCBI
#dt = table(H3_NCBI$Year, H3_NCBI$BirthYear); use.cols = which(as.character(1990:2014) %in% rownames(dt))
bys = min(H3_NCBI$BirthYear):2010; yrs = sort(unique(H3_NCBI$Year))
dt = matrix(NA, length(yrs), length(bys), dimnames = list(yrs, bys))
for(yy in 1:length(yrs)){
  for(ii in 1:length(bys)){
    dt[yy, ii] = sum(H3_NCBI$BirthYear == bys[ii] & H3_NCBI$Year == yrs[yy])
  }}
use.cols = which(as.character(1990:2014) %in% rownames(dt))
xx = barplot(dt, col = cols[use.cols], border = NA, main = 'H3 NCBI', space = 0, xaxt = 'n')
use.ax = which(colnames(dt) %in% as.character(yys))
nms = colnames(dt)[use.ax]
axis(1, at = xx[use.ax], labels = nms)




plot.new(); par(mfrow = c(1,1))
plot.window(xlim = c(0,1), ylim = c(0,1))
legend('left', legend = as.character(1990:2002), fill = cols[1:13])
legend('right', legend = as.character(2003:2014), fill = cols[14:24])






### Plot data vs. age
## H1 AZ
par(mfrow = c(2,2))
dt = table(H1_AZ$Year, H1_AZ$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H1 Arizona')

## H1 NCBI
dt = table(H1_NCBI$Year, H1_NCBI$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H1 NCBI')

## H3 AZ
dt = table(H3_AZ$Year, H3_AZ$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H3 Arizona')

## H3 NCBI
dt = table(H3_NCBI$Year, H3_NCBI$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H3 NCBI')

table(H1_AZ$year)
table(H1_NCBI$Year)



## Exclude seasons with less than 100 cases
threshold = 100
c(sum(table(H1_AZ$SEASON) > threshold), 
sum(table(H3_AZ$SEASON) > threshold),
sum(table(H1_NCBI$Year) > threshold),
sum(table(H3_NCBI$Year) > threshold))



## Trim off years that don't have enough cases to meet the threshold
## Other cases can be used in the test set later
keep = names(which(table(H1_AZ$SEASON) > threshold))
H1_AZ_trimmed = H1_AZ[H1_AZ$SEASON %in% keep, ]

keep = names(which(table(H3_AZ$SEASON) > threshold))
H3_AZ_trimmed = H3_AZ[H3_AZ$SEASON %in% keep, ]

keep = names(which(table(H1_NCBI$Year) > threshold))
H1_NCBI_trimmed = H1_NCBI[H1_NCBI$Year %in% keep, ]

keep = names(which(table(H3_NCBI$Year) > threshold))
H3_NCBI_trimmed = H3_NCBI[H3_NCBI$Year %in% keep, ]


## Write a function to subsample 100 observations from each year of observation
# balanced.data = function(NCBI.dat.in, threshold){
#   years = unique(NCBI.dat.in$Year)
#   valid = which(NCBI.dat.in$Year == years[1])
#   subsample = NCBI.dat.in[sample(valid, size = threshold, replace = FALSE),]
#   for(ii in 1:length(years)){
#     valid = which(NCBI.dat.in$Year == years[ii])
#     subsample = rbind(subsample, 
#                   NCBI.dat.in[sample(valid, size = threshold, replace = FALSE),])
#   }
#   subset(subsample, select = c("Length", "Host", "Protein", "Subtype", "Country", "Region", "Year", "Virus", "Age", "BirthYear"))
# }


## Write a function to bootstrap cases from the combined raw data
##    Sample each point with a probability inversely proportional to the year's representation in the data
##    This balances the data with respect to year sampled.
##    Sample with repeition the same number of cases observed in H1 and H3 data
##    Combine the H1 and H3 samples into one big data set
boot.age = function(){
  nn.H1 = nrow(H1_NCBI_trimmed)
  nn.H3 = nrow(H3_NCBI_trimmed)
  
  ## Weight sampling prob by number obs/year
  weights.H1 = 1/table(H1_NCBI_trimmed$Year)
  weights.H3 = 1/table(H3_NCBI_trimmed$Year)
  
  ## Create weights vector
  weights.vec.H1 = weights.H1[as.character(H1_NCBI_trimmed$Year)]
  weights.vec.H3 = weights.H3[as.character(H3_NCBI_trimmed$Year)]
  
  ## Bootstrap, sampling each point with a probability inversely proportional to that year's representation in the dataset
  boot.ind.H1 = sample(x = 1:nn.H1, size = nn.H1, replace = TRUE, prob = weights.vec.H1)
  boot.ind.H3 = sample(x = 1:nn.H3, size = nn.H1, replace = TRUE, prob = weights.vec.H3) ## Sample nn.H1 of each so that there are an equal # of H1 and H3 cases represented
  return(rbind(H1_NCBI_trimmed[boot.ind.H1, ], H3_NCBI_trimmed[boot.ind.H3, ]))
}


# Bootstrap null
boot.null = boot.age()
# Calculate the fraciton of each age group in the bootstrapped null
null.ages = table(boot.null$Age)/nrow(boot.null); null.ages

## Plot the bootstrapped null vs. observed age dist
par(mfrow = c(1,1))
plot(as.numeric(names(null.ages)), null.ages, pch = 16)
H1.ages = table(H1_AZ_trimmed$Age)/nrow(H1_AZ_trimmed)
points(as.numeric(names(H1.ages)), H1.ages, col = 'blue')
H3.ages = table(H3_AZ_trimmed$Age)/nrow(H3_AZ_trimmed)
points(as.numeric(names(H3.ages)), H3.ages, col = 'red')
  
## Calculate a smoothing spline for the null values
library(splines)
age.spline = smooth.spline(null.ages, keep.data = TRUE,  df = 10)
spline.vals = predict(age.spline, x = 0:120)
spline.vals$y[which(spline.vals$y < 0)] = 0
spline.vals$y[1] = .05
par(mfrow = c(1,1))
plot(y~x, data = spline.vals, type = 'b', ylim = c(0, .06))
# Compare to actual values in the bootstrapped data.
points(as.numeric(names(null.ages)), null.ages, col = 'green', pch = names(null.ages))
points(as.numeric(names(H1.ages)), H1.ages, col = 'blue')
points(as.numeric(names(H3.ages)), H3.ages, col = 'red')



################# 
# Write a funciton to reformat data into a matrix with birth years 1918:2017 on the columns and the year of observation on the rows
reformat.data = function(dat.in){
  yrs = unique(dat.in$Year)
  bys = 1918:2017
  out.mat = matrix(NA, nrow = length(yrs), ncol = length(bys), dimnames = list(yrs, bys))
  for(yy in 1:length(yrs)){
    for(bb in 1:length(bys)){
      out.mat[yy, bb] = sum(dat.in$Year == yrs[yy] & dat.in$BirthYear == bys[bb])
    }
  }
  out.mat
}



## Reformatted data inputs
H1.inputs = reformat.data(H1_AZ_trimmed)
H3.inputs = reformat.data(H3_AZ_trimmed)



## Age inputs
## Reformat age inputs to the same dimensions as the H1 and H3 data
age.in.H1 = H1.inputs*0
for(ii in 1:nrow(age.in.H1)){
  yr = as.numeric(rownames(age.in.H1)[ii]) # Set the yr of observation
  n.valid.ages = min(yr-1918+1, ncol(age.in.H1)) # From the perspective of yr, how many age groups were born since 1918?
  age.in.H1[ii, 1:n.valid.ages] = (spline.vals$y[n.valid.ages:1]) # Fill in the desired ages with null spline values
}

age.in.H3 = H3.inputs*0
for(ii in 1:nrow(age.in.H3)){
  yr = as.numeric(rownames(age.in.H3)[ii])
  n.valid.ages = min(yr-1918+1, ncol(age.in.H3))
  age.in.H3[ii, 1:n.valid.ages] = (spline.vals$y[n.valid.ages:1])
}


## Assume risk decreases with age? Come back to this.
# exp.age.decrease = exp(seq(0, log(.4), length = 120))
# linear.age.decrease = seq(.4, 1, length = 120)
# 
# linear.age.decrease.H3 = exp.age.decrease.H3 = H3.inputs*0
# for(ii in 1:nrow(age.in.H3)){
#   yr = as.numeric(rownames(age.in.H3)[ii])
#   n.valid.ages = min(yr-1918+1, ncol(age.in.H3))
#   linear.age.decrease.H3[ii, 1:n.valid.ages] = (linear.age.decrease[1:n.valid.ages])
#   exp.age.decrease.H3[ii, 1:n.valid.ages] = (linear.age.decrease[1:n.valid.ages])
# }
# 
# linear.age.decrease.H1 = exp.age.decrease.H1 = H1.inputs*0
# for(ii in 1:nrow(age.in.H1)){
#   yr = as.numeric(rownames(age.in.H1)[ii])
#   n.valid.ages = min(yr-1918+1, ncol(age.in.H1))
#   linear.age.decrease.H1[ii, 1:n.valid.ages] = (linear.age.decrease[1:n.valid.ages])
#   exp.age.decrease.H1[ii, 1:n.valid.ages] = (linear.age.decrease[1:n.valid.ages])
# }


## Reconstruct imprinting patterns
source('Infection.age.structure.vaccination.R') ## This script reconstructs patterns for USA and other countires up to 2017, assmes vaccination delays imprinting
vaccination.matrix.raw = read.csv('~/Dropbox/R/ImmuneAgeStructure/ForSimulations/Demography_files/vaccinationRates_ages0to9.csv')
vaccination.matrix = cbind(year = vaccination.matrix.raw$year, vaccination.matrix.raw[,-1]*0.6) #Assume 0.6 in naiive children


## Chop off older years with no data
H1.min = which.max(colSums(H1.inputs) > 0); H1.inputs = H1.inputs[, -(1:(H1.min-1))]
## Chop off the oldest case for H1, since it's way older than others
H1.min = which.max(colSums(H1.inputs) > 0); H1.inputs = H1.inputs[, -(1:(H1.min-1))]

H3.min = which.max(colSums(H3.inputs) > 0); H3.inputs = H3.inputs[, -(1:(H3.min-1))]

max.by = max(c(as.numeric(rownames(H1.inputs)), as.numeric(rownames(H3.inputs))))
H1.inputs = H1.inputs[, -which(as.numeric(colnames(H1.inputs)) > max.by)]
H3.inputs = H3.inputs[, -which(as.numeric(colnames(H3.inputs)) > max.by)]

## Get row and column years so as to arrange weights to match data
use.years.H1 = gsub("(\\d)(\\d)(\\d)(\\d)(\\w+)", "\\1", rownames(H1.inputs))
use.birth.years.H1 = gsub("(\\d)(\\d)(\\d)(\\d)(\\w+)", "\\1", colnames(H1.inputs))
use.years.H3 = gsub("(\\d)(\\d)(\\d)(\\d)(\\w+)", "\\1", rownames(H3.inputs))
use.birth.years.H3 =gsub("(\\d)(\\d)(\\d)(\\d)(\\w+)", "\\1", colnames(H3.inputs))
use.years = unique(c(use.years.H1, use.years.H3))

wts = get.type.weights.AB.vaccination_2(years.out = as.numeric(use.years), Countries.out = c('USA'), type = 5, vax.matrix = vaccination.matrix)

H1.rows = paste(use.years.H1, 'USA', sep = '')
H3.rows = paste(use.years.H3, 'USA', sep = '')

wm.H1 = wts$weights.master.1[H1.rows, use.birth.years.H1] + wts$weights.master.2[H1.rows, use.birth.years.H1]
wo.H1 = wts$weights.master.3[H1.rows, use.birth.years.H1] + wts$weights.master.naiive[H1.rows, use.birth.years.H1]

wm.H1.subtype = wts$weights.master.1[H1.rows, use.birth.years.H1] 
wo.H1.subtype = wts$weights.master.3[H1.rows, use.birth.years.H1] + wts$weights.master.naiive[H1.rows, use.birth.years.H1] + wts$weights.master.2[H1.rows, use.birth.years.H1]

wm.H3 = wts$weights.master.3[H3.rows, use.birth.years.H3]
wo.H3 = wts$weights.master.1[H3.rows, use.birth.years.H3] + wts$weights.master.2[H3.rows, use.birth.years.H3]+ wts$weights.master.naiive[H3.rows, use.birth.years.H3]

## Weights that include NA patterns
H1N1.wts.H1 = wts$weights.master.1[H1.rows, use.birth.years.H1]
H2N2.wts.H1 = wts$weights.master.2[H1.rows, use.birth.years.H1]
H3N2.wts.H1 = wts$weights.master.3[H1.rows, use.birth.years.H1]
naive.wts.H1 = wts$weights.master.naiive[H1.rows, use.birth.years.H1]

H1N1.wts.H3 = wts$weights.master.1[H3.rows, use.birth.years.H3]
H2N2.wts.H3 = wts$weights.master.2[H3.rows, use.birth.years.H3]
H3N2.wts.H3 = wts$weights.master.3[H3.rows, use.birth.years.H3]
naive.wts.H3 = wts$weights.master.naiive[H3.rows, use.birth.years.H3]


## ADD NA PATTERNS HERE AND BELOW!

setwd('~/Dropbox/R/2017_seasonal_flu/')

age.in.H1 = age.in.H1[use.years.H1, use.birth.years.H1]
age.in.H3 = age.in.H3[use.years.H3, use.birth.years.H3]
# linear.age.decrease.H1 = linear.age.decrease.H1[use.years.H1, use.birth.years.H1]
# linear.age.decrease.H3 = linear.age.decrease.H3[use.years.H3, use.birth.years.H3]




## Fit models!
source('Seasonal_likelihood_functions.R')

## No age effects
age.mod.H1 = lk.A.H(pars = c(Hm = 1), wm = wm.H1, wo = wo.H1, age.spline = age.in.H1, xx = H1.inputs); age.mod.H1
# Age and group-level imprinting
age.by.mod.H1 = optim(par = c(Hm = .2), fn = lk.A.H,  wm = wm.H1, wo = wo.H1, age.spline = age.in.H1, xx = H1.inputs, method = 'L-BFGS-B', lower = 0.001, upper = 1); age.by.mod.H1
# Age and subtype-level imprinting
age.by.subtype.mod.H1 = optim(par = c(Hm = .2), fn = lk.A.H,  wm = wm.H1.subtype, wo = wo.H1.subtype, age.spline = age.in.H1, xx = H1.inputs, method = 'L-BFGS-B', lower = 0.001, upper = 1); age.by.subtype.mod.H1
# Age and NA imprinting
age.NA.mod.H1 = optim(par = c(Nm = .2), fn = lk.A.N,  wm = H1N1.wts.H1, wo = H2N2.wts.H1+H3N2.wts.H1+naive.wts.H1, age.spline = age.in.H1, xx = H1.inputs, method = 'L-BFGS-B', lower = 0.001, upper = 1); age.NA.mod.H1
# Age and NA and HA imprinting
age.HA.NA.mod.H1 = optim(par = c(Hm = .5, Nm = .5), fn = lk.A.N.H, wmm = H1N1.wts.H1, wmo = H2N2.wts.H1, wom = H2N2.wts.H1*0, woo = H3N2.wts.H1+naive.wts.H1, age.spline = age.in.H1, xx = H1.inputs, method = 'L-BFGS-B', lower = .001, upper = 1); age.HA.NA.mod.H1

# Age and group-level imprinting
#age.by.int.mod.H1 = optim(par = c(Hm = .2), fn = lk.A.H.I,  wm = wm.H1, wo = wo.H1, age.spline = age.in.H1, age.decrease = linear.age.decrease.H1, xx = H1.inputs, method = 'L-BFGS-B', lower = 0.001, upper = 1); age.by.int.mod.H1
# Age and subtype-level imprinting
#age.by.subtype.int.mod.H1 = optim(par = c(Hm = .2), fn = lk.A.H.I,  wm = wm.H1.subtype, wo = wo.H1.subtype, age.decrease = linear.age.decrease.H1, age.spline = age.in.H1, xx = H1.inputs, method = 'L-BFGS-B', lower = 0.001, upper = 1); age.by.subtype.int.mod.H1


## Calculate AIC
AIC = function(nll, n.pars){
  2*n.pars+2*nll
}

AIC.H1.NCBI = AIC(c(age.mod.H1, age.by.mod.H1$value, age.by.subtype.mod.H1$value, age.NA.mod.H1$value, 
                    age.HA.NA.mod.H1$value), c(0, 1, 1, 1, 2)); AIC.H1.NCBI
names(AIC.H1.NCBI) = c('A', 'AH', 'AHs', 'AN', 'ANH')
# del.aic
del.AIC.H1.NCBI = sort(AIC.H1.NCBI-min(AIC.H1.NCBI)); del.AIC.H1.NCBI

# Calculate and scale predictions
A.prediction = age.in.H1/rowSums(age.in.H1)*rowSums(H1.inputs)
# AH
p.raw.AH = ( age.in.H1*(age.by.mod.H1$par['Hm']*wm.H1+wo.H1)  )
AH.prediction = p.raw.AH/rowSums(p.raw.AH)*rowSums(H1.inputs)
# AN = AHs
p.raw.AHs= ( age.in.H1*(age.by.subtype.mod.H1$par['Hm']*wm.H1.subtype+wo.H1.subtype)  )
AHs.prediction = p.raw.AHs/rowSums(p.raw.AHs)*rowSums(H1.inputs)
# ANH
hh = age.HA.NA.mod.H1$par['Hm']
nn = age.HA.NA.mod.H1$par['Nm']
p.raw.ANH= ( age.in.H1*(hh*nn*H1N1.wts.H1 + hh*H2N2.wts.H1 + H3N2.wts.H1+naive.wts.H1)  )
ANH.prediction = p.raw.ANH/rowSums(p.raw.ANH)*rowSums(H1.inputs)




## Plot results
par(mfrow = c(1,2))
xx = as.numeric(use.birth.years.H1)
xx = barplot(colSums(H1.inputs), main = 'H1N1, NCBI Data', xlab = 'birth year', ylab = 'case count', col = 'gray', border = NA, space = 0, xaxt = 'n')
axis(1, at = which(colnames(H1.inputs) %in% seq(1930, 2000, by = 10)), labels = colnames(H1.inputs)[which(colnames(H1.inputs) %in% seq(1930, 2000, by = 10))])
lines(xx, colSums(A.prediction), col = 'red', lwd = 2)
lines(xx, colSums(AH.prediction), col = 'blue', lwd = 2)
lines(xx, colSums(AHs.prediction), col = 'navy', lwd = 2)
lines(xx, colSums(ANH.prediction), col = 'yellow', lwd = 2)
legend('topleft', c('Age only', 'Age + group-level imprinting', 'Age + subtype-level imprinting'), lty = 1, col = c('red', 'blue', 'navy'), bty = 'n')




### REPEAT FOR H3
## No age effects
age.mod.H3 = lk.A.H(pars = c(Hm = 1), wm = wm.H3, wo = wo.H3, age.spline = age.in.H3, xx = H3.inputs); age.mod.H3
# Age plus birth year
age.by.mod.H3 = optim(par = c(Hm = .5), fn = lk.A.H,  wm = wm.H3, wo = wo.H3, age.spline = age.in.H3, xx = H3.inputs, method = 'L-BFGS-B', lower = .001, upper = 1); age.by.mod.H3
# Age, NA imprinting
age.NA.mod.H3 = optim(par = c(Nm = .2), fn = lk.A.N,  wm = H2N2.wts.H3+H3N2.wts.H3, wo = H1N1.wts.H3+naive.wts.H3, age.spline = age.in.H3, xx = H3.inputs, method = 'L-BFGS-B', lower = .001, upper = 1); age.NA.mod.H3
# Age, HA, NA imprinting
age.NA.HA.mod.H3 = optim(par = c(Hm = .9, Nm = .9), fn = lk.A.N.H,  wmm = H3N2.wts.H3, wmo = H1N1.wts.H3*0, wom = H2N2.wts.H3, woo = naive.wts.H3+H1N1.wts.H3, age.spline = age.in.H3, xx = H3.inputs, method = 'L-BFGS-B', lower = .001, upper = 1); age.NA.HA.mod.H3

# Age plus birth year, decreasing with age
# Age plus birth year
# age.by.int.mod.H3 = optim(par = c(Hm = .2), fn = lk.A.H.I,  wm = wm.H3, wo = wo.H3, age.spline = age.in.H3, age.decrease = linear.age.decrease.H3, xx = H3.inputs); age.by.int.mod.H3



AIC.H3.NCBI = AIC(c(age.mod.H3, age.by.mod.H3$value, age.NA.mod.H3$value, age.NA.HA.mod.H3$value), c(0, 1, 1, 2)); AIC.H3.NCBI
names(AIC.H3.NCBI) = c('A', 'AH', 'AN', 'ANH')
# del.aic
del.AIC.H3.NCBI = AIC.H3.NCBI-min(AIC.H3.NCBI); sort(del.AIC.H3.NCBI)



## Overall AIC based on sum of AIC from fits to H1 and H3
AIC.combined = AIC(c(age.mod.H3+age.mod.H1, 
                     age.by.mod.H3$value+age.by.mod.H1$value, 
                     age.NA.mod.H3$value+age.NA.mod.H1$value), c(0, 2, 2)); AIC.combined
names(AIC.combined) = c('A', 'AH', 'AN')
# del.aic
del.AIC.combined = AIC.combined-min(AIC.combined); sort(del.AIC.combined)




# Calculate and scale predictions
# Calculate and scale predictions
A.prediction = age.in.H3/rowSums(age.in.H3)*rowSums(H3.inputs)
# AH
p.raw.AH = ( age.in.H3*(age.by.mod.H3$par['Hm']*wm.H3+wo.H3)  )
AH.prediction = p.raw.AH/rowSums(p.raw.AH)*rowSums(H3.inputs)
# AN 
p.raw.AN= ( age.in.H3*(age.NA.mod.H3$par['Nm']*H3N2.wts.H3 + H1N1.wts.H3+H2N2.wts.H3+naive.wts.H3)  )
AN.prediction = p.raw.AN/rowSums(p.raw.AN)*rowSums(H3.inputs)
# ANH
hh = age.NA.HA.mod.H3$par['Hm']
nn = age.NA.HA.mod.H3$par['Nm']
p.raw.ANH= ( age.in.H3*(hh*nn*H3N2.wts.H3 + nn*H2N2.wts.H3 + H1N1.wts.H3+naive.wts.H3)  )
ANH.prediction = p.raw.ANH/rowSums(p.raw.ANH)*rowSums(H3.inputs)



## Plot results
xx = as.numeric(use.birth.years.H3)
xx = barplot(colSums(H3.inputs), main = 'H3N2, NCBI Data', xlab = 'birth year', ylab = 'case count', col = 'gray', border = NA, xaxt = 'n', space = 0)
axis(1, at = which(colnames(H3.inputs) %in% seq(1930, 2000, by = 10)), labels = colnames(H3.inputs)[which(colnames(H3.inputs) %in% seq(1930, 2000, by = 10))])
lines(xx, colSums(A.prediction), col = 'red', lwd = 2)
lines(xx, colSums(AH.prediction), col = 'blue', lwd = 2)
lines(xx, colSums(AN.prediction), col = 'navy', lwd = 2)
lines(xx, colSums(ANH.prediction), col = '#e1ff00', lwd = 2)
legend('topleft', c('Age only', 'Age + group-level imprinting', 'Age + subtype-level imprinting'), lty = 1, col = c('red', 'blue', 'navy'), bty = 'n')






del.AIC.H1.NCBI
sort(del.AIC.H3.NCBI)







#### Plot age and birth year distributions over time for H3N2
H3N2.dat = data.frame(rbind(cbind(H3_NCBI$Year, H3_NCBI$BirthYear, H3_NCBI$Age),
                     cbind(H3_AZ$Year, H3_AZ$BirthYear, H3_AZ$Age)))
colnames(H3N2.dat) = c('Year', 'BirthYear', 'Age')

table(H3N2.dat$Year)

pdf('Patterns_over_time_H3N2.pdf', height = 16)
n.years = length(table(H3N2.dat$Year))
#years = sort(unique(H3N2.dat$Year))
years = c(1993, 2000, 2005, 2010, 2014)
n.years = length(years)
par(mfrow = c(n.years, 2))
for(ii in 1:n.years){
  sub = subset(H3N2.dat, Year == years[ii])
  hist(sub$Age, breaks = c(seq(0, 100, by = 5), 150), main = paste(years[ii], '   N = ', nrow(sub), sep = ''), xlab = 'Age', ylab = '', freq = TRUE)
  abline(v = years[ii] - c(1968, 1957), col = 'red')
  hist(sub$BirthYear, breaks = c(1887, 1918, 1957, 1968, 1977, 2009, 2017), main = paste(years[ii], '   N = ', nrow(sub), sep = ''), xlab = 'Birth Year', ylab = '', freq = TRUE)
}
dev.off()
  



H1N1.dat = data.frame(rbind(cbind(H1_NCBI$Year, H1_NCBI$BirthYear, H1_NCBI$Age),
                            cbind(H1_AZ$Year, H1_AZ$BirthYear, H1_AZ$Age)))
colnames(H1N1.dat) = c('Year', 'BirthYear', 'Age')

table(H1N1.dat$Year)

pdf('Patterns_over_time_H1N1.pdf', height = 12)
#years = sort(unique(H1N1.dat$Year))
years = c(2001, 2007, 2010, 2013)
n.years = length(years)
par(mfrow = c(n.years, 2))
for(ii in 1:n.years){
  sub = subset(H1N1.dat, Year == years[ii])
  hist(sub$Age, breaks = c(seq(0, 100, by = 5), 150), main = paste(years[ii], '   N = ', nrow(sub), sep = ''), xlab = 'Age', ylab = '', freq = TRUE)
  abline(v = years[ii] - c(1968, 1957), col = 'red')
  hist(sub$BirthYear, breaks = c(1887, 1918, 1957, 1968, 1977, 2009, 2017), main = paste(years[ii], '   N = ', nrow(sub), sep = ''), xlab = 'Birth Year', ylab = '', freq = TRUE)
}
dev.off()



old.H3 = sum(H3N2.dat$Age >= 65)
old.H1 = sum(H1N1.dat$Age >= 65)

old.H3/old.H1


old.H3 = sum(H3N2.dat$BirthYear <= 1957)
old.H1 = sum(H1N1.dat$BirthYear <= 1957)

old.H3/old.H1
