#### Visualize two kinds of sequence data and test GAMMS
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')

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


## Calculate age
names(raw.dat) = c('season', 'birthyear', 'subtype')
## Exclude cases born before 1918 (imprinting history not known, and there are few cases)
raw.dat = raw.dat[-which(raw.dat$birthyear < 1918), ]
raw.dat$year1 = as.numeric(gsub(pattern = '(\\d\\d\\d\\d)\\d\\d', replacement = "\\1", x = raw.dat$season))
raw.dat$year2 = raw.dat$year1 + 1
raw.dat$age = raw.dat$year2 - raw.dat$birthyear



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
# source('Infection.age.structure.vaccination.R') ## This script reconstructs patterns for USA and other countires up to 2017, assmes vaccination delays imprinting
# vaccination.matrix.raw = read.csv('~/Dropbox/R/ImmuneAgeStructure/ForSimulations/Demography_files/vaccinationRates_ages0to9.csv')
# vaccination.matrix = cbind(year = vaccination.matrix.raw$year, vaccination.matrix.raw[,-1]*0.6) #Assume 0.6 in naiive children
# wts = get.type.weights.AB.vaccination_2(years.out = as.numeric(use.years), Countries.out = c('USA'), type = 5, vax.matrix = vaccination.matrix)


# source('../Reconstructions/Infection_age_structure_2017.R')
# weights = get.type.weights.AB(years.out = c(1994, 1995, 2003:2015), Countries.out = 'USA', region.in = c('default'))
# save(weights, file = 'AZ_weights.RData')
load('AZ_weights.RData')

proH1.master = weights[[1]][,as.character(2015:1918)]
proH3.master = weights[[3]][,as.character(2015:1918)]
prog1.master = weights[[1]][,as.character(2015:1918)]+weights[[2]][,as.character(2015:1918)]
prog2.master = proH3.master
proN1.master = proH1.master
proN2.master = weights[[2]][,as.character(2015:1918)]+weights[[3]][,as.character(2015:1918)]




#### Define age-specific indicators
agemat = t(sapply(X = c(1994, 1995, 2003:2015), function(xx) xx - 2015:1918))
rownames(agemat) = ssns; colnames(agemat) = bys
inclusion = matrix(1, nrow(agemat), ncol(agemat))
## Exclude entries over age 73 from predicted probabilities and data.
inclusion[which(agemat>65)] = 0
H1.master = H1.master*inclusion
H3.master = H3.master*inclusion

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

barplot(a18.24)
barplot(a25.31)
barplot(a32.38)
barplot(a39.45)
barplot(a46.52)
barplot(a53.59)
barplot(a60.66)
barplot(a67.73)
barplot(a74.80)
barplot(a81.90plus)
barplot(a0.4 + a5.10 + a11.17 + a18.24+a25.31+a32.38+a39.45+a46.52+a53.59+a60.66+a67.73+a74.80+a81.90plus, col = rainbow(14))





#######################################
## Define likelihood
######################################
#### INPUTS:
####    pars - named vector of pars to be fit. If not all four pars exist in a given model, code automatically assigns the null value of 1. Parameters not named in "pars" will not be fit.
####    wPro.H1 - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_i) who are protected against H1 based on childhood imprinting
####    wPro.H3 - same as above, but characterizes protection against seasonal H3N2
####    dat.H1 - n-vector or mxn matrix of observed case counts. Each entry represents the number of cases of a given age (n_i), observed in a given country and season (m_i).
####    a18.24- n-vector or mxn matrix of indicator variables. 1 indicates membership in the 18-24 age class, 0 otherwise
####    All other a##.## inputs are indicators taht follow the same logic.

##### OUTPUTS:
#####    negative log likelihood (scalar)


##### NOTES:  
####    All models tested are nested in this likelihood structure. E.g. to exclude the effects of protection from the model, input wPro = 0, and omit "rPro" from pars.

nll = function(pars, wPro.H1, dat.H1, wPro.H3, dat.H3, a0.4, a5.10, a11.17, a18.24, a25.31, a32.38, a39.45, a46.52, a53.59, a60.66, a67.73, a74.80, a81.90plus){
  # 1. Assign parameters to be fit
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  b = 1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r5.10 = pars['r5.10'] # Expected risk for 5 to 10 year olds
  r11.17 = pars['r11.17'] # etc.
  r18.24 = pars['r18.24'] 
  r25.31 = pars['r25.31']
  r32.38 = pars['r32.38']
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73'] 
  r74.80 = pars['r74.80'] 
  r81.90p = pars['r81.90p'] 
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.baseline = (b*a0.4 +b*r5.10*a5.10 +b*r11.17*a11.17+b*r18.24*a18.24+b*r25.31*a25.31+ b*r32.38*a32.38+ b*r39.45*a39.45+ b*r46.52*a46.52+ b*r53.59*a53.59+ b*r60.66*a60.66+ b*r67.73*a67.73+ b*r74.80*a74.80+ b*r81.90p*a81.90plus)
  age.baseline = age.baseline/rowSums(age.baseline)
  
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = inclusion*age.baseline * (wPro.H1*rPro.H1+(1-wPro.H1))
  
  pp.H3 = inclusion*age.baseline * (wPro.H3*rPro.H3+(1-wPro.H3))
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(dat.H1))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
    lk.H1 = -dmultinom(dat.H1, size = sum(dat.H1), prob = pp.H1, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(dat.H1)[1])
    for(jj in 1:dim(dat.H1)[1]){ #Find the neg log density for each row (dim 1) and take the sum
      storage[jj] = -dmultinom(dat.H1[jj,], size = sum(dat.H1[jj,]), prob = pp.H1[jj,], log = TRUE)
    }
    lk.H1 = sum(storage) 
  }
  if(is.null(dim(dat.H3))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
    lk.H3 = -dmultinom(dat.H3, size = sum(dat.H3), prob = pp.H3, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(dat.H3)[1])
    for(jj in 1:dim(dat.H3)[1]){ #Find the neg log density for each row (dim 1) and take the sum
      storage[jj] = -dmultinom(dat.H3[jj,], size = sum(dat.H3[jj,]), prob = pp.H3[jj,], log = TRUE)
    }
    lk.H3 = sum(storage) 
  }
  lk.H1+lk.H3 # end function 
}

# Test function
nll(pars = c(rPro.H1 = .5, rPro.H3 = .5, r5.10 = 1.1, r11.17 = .9, r18.24 = .9, r25.31 = .9, r32.38 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73= .9, r74.80 = .9, r81.90p = .9), wPro.H1 = proH1.master, wPro.H3 = proH3.master, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, dat.H1 = H1.master, dat.H3 = H3.master)






## Write a function wrapper, so that you only have to input a vector of initial par values and par value limits involved in model comparison
nll.wrapper = function(pars.in, pro.H1, pro.H3, lower.in, upper.in){
  ## Concatenate vector to initialize pars that may or may not be included, with a vector of all age pars (always included)
  pvec = c(pars.in, r5.10 = 1.1, r11.17 = .9, r18.24 = .9, r25.31 = .9, r32.38 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73= .9, r74.80 = .9, r81.90p = .9)
  
  optim(par = pvec, fn = nll, wPro.H1 = pro.H1, wPro.H3 = pro.H3, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, dat.H1 = H1.master, dat.H3 = H3.master, method = 'L-BFGS-B', lower = c(lower.in, rep(.001, 12)), upper = c(upper.in, rep(5, 12)))
}

## Test
## Maximal model
nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(.001, .001), upper.in = c(1, 1))

## Reduced model, vaccination and imprinting only
nll.wrapper(pars.in = NULL, pro.H1 = 1, pro.H3 = 1, lower.in = NULL, upper.in = NULL)



#######################################
## Model comparison 
##      null is all cases tested
##      fit to H1N1
######################################
pro.low = .001; pro.high = 1 # Upper and lower bounds for protection estimates
## set upper and lower bounds for all pars
## 29. G 
lk.G = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.G

## 30. S
lk.S = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.S

## 31. N
lk.N = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.N

## 32. NULL
lk.null = nll.wrapper(pars.in = NULL, pro.H1 = 0, pro.H3 = 0, lower.in = NULL, upper = NULL); lk.null






## Pull out the variable names that store likelihoods
mods = mget(ls(pattern = "lk."))
nll = numeric(length(mods))
AICs = numeric(length(mods))
for(ii in 1:length(mods)){
  nll[ii] = mods[[ii]]$value
  AICs[ii] = 2*length(mods[[ii]]$par)+2*mods[[ii]]$value
}

names(AICs) = names(mods)
names(nll) = names(mods)
AICs = sort(AICs)
del.AIC = AICs - min(AICs)
del.AIC



## Remove degenerate models
ii = 1
while(ii < length(del.AIC)){
  remainder = del.AIC[(ii+1):length(del.AIC)] - del.AIC[ii]
  drop = which(round(remainder, 3) %% 2 == 0)
  del.AIC = del.AIC[!names(del.AIC) %in% names(drop)]
  ii = ii+1
}
del.AIC







## Plot model results
plotmod1 = function(pars, pro.H1 = 1, pro.H3 = 1, i.type = NULL){
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  b = 1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r5.10 = pars['r5.10'] # Expected risk for 5 to 10 year olds
  r11.17 = pars['r11.17'] # etc.
  r18.24 = pars['r18.24'] 
  r25.31 = pars['r25.31']
  r32.38 = pars['r32.38']
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73'] 
  r74.80 = pars['r74.80'] 
  r81.90p = pars['r81.90p'] 
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.baseline = (b*a0.4 +b*r5.10*a5.10 +b*r11.17*a11.17+b*r18.24*a18.24+b*r25.31*a25.31+ b*r32.38*a32.38+ b*r39.45*a39.45+ b*r46.52*a46.52+ b*r53.59*a53.59+ b*r60.66*a60.66+ b*r67.73*a67.73+ b*r74.80*a74.80+ b*r81.90p*a81.90plus)
  age.baseline = age.baseline/rowSums(age.baseline)
  
  imprinting.H1 = (pro.H1*rPro.H1+(1-pro.H1))
  imprinting.H3 = (pro.H3*rPro.H3+(1-pro.H3))
  
  par(mfrow = c(2,2))
  ## All rows in age baseline are the same, so plot one arbitrarily
  plot(0:97, age.baseline[15,], main = 'Age effects', ylab = 'proportion of cases', xlab = 'case age')
  
  
  ## Imprinting protection
  if(is.na(pars['rPro.H1'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste('Imprinting protection\n', i.type, sep = ''))
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(colMeans(imprinting.H1), col = 'dodgerblue', cex = .7, ylab = 'est. relative risk impact', xlab = 'case birth year', main = 'Imprinting', ylim = c(0,1), xaxt = 'n')
    points(colMeans(imprinting.H3), col = 'firebrick1', cex = .7)
    axis(1, at = seq(3, 98, by = 10), labels = seq(1920, 2015, by = 10), las = 2) }
  
  
  ## Relative risk
  parnames = c("ages 0-4", "ages5-10", "ages 11-17", "ages 18-24",  "ages 25-31", "ages 32-38",  "ages 39-45",  "ages 46-52",  "ages 53-59",  "ages 60-66", "ages 67-73",  "ages 74-80",  "ages 81+", "", "", "", paste('imprinting', i.type, sep = ', ')) 
  xvals = c(pars[grep(pattern = "r\\d\\d?.\\d\\d", x = names(pars))], NA, NA, NA, pars['rPro.H1'], pars['rPro.H3'], 1)
  yvals = c(1:15, 16, 16, 0)+1
  par(mar = c(4, 7, 2, 1))
  plot(xvals, yvals, xlim = c(0, 1.3), xaxt = 'n', yaxt = 'n', xlab = 'Relative risk estimate', ylab = '', col = c(rep('black', 12), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), pch = 9)
  axis(side = 1, at = seq(0, 1.25, by = .25))
  axis(side = 2, at = 1:17, labels = parnames, las = 2)
  abline(v = 1, lty = 2)
  
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = inclusion*age.baseline * imprinting.H1
  pp.H3 = inclusion*age.baseline * imprinting.H3
  return(rbind(colSums(pp.H1/rowSums(pp.H1)*rowSums(H1.master)), colSums(pp.H3/rowSums(pp.H3)*rowSums(H3.master))))
}


## Plot best model
npred = plotmod1(lk.N$par, pro.H1 = proN1.master, pro.H3 = proN2.master, i.type = 'NA Subtype')
spred = plotmod1(lk.S$par, pro.H1 = proH1.master, pro.H3 = proH3.master, i.type = 'HA Subtype')

pdf('AZ_predictions.pdf', width = 7, height = 4.5)
par(mfrow = c(1,2))
plot(2015:1918, colSums(H1.master), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H1N1')
lines(2015:1918, npred[1,], col = 'magenta', lwd = 2)
lines(2015:1918, spred[1,], col = 'orange', lwd = 2)
legend('topleft', legend = c('N', 'S'), col = c('magenta', 'orange'), lty = 1)

plot(2015:1918, colSums(H3.master), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H3N2')
lines(2015:1918, npred[2,], col = 'magenta', lwd = 2)
lines(2015:1918, spred[2,], col = 'orange', lwd = 2)
dev.off()


rr = function(center, width = 1, height = 1, col.in = 'navy'){
  rect(xleft = center[1]-width/2, ybottom = center[2]-height/2, xright = center[1]+width/2, ytop = center[2]+height/2, col = col.in, border = 'black')
}

dev.off()
pdf('AZ_AIC.pdf', width = 4, height = 3)
## Repeat, only plotting non-degenerate models
## Plot AIC results
mods2 = mods[rev(names(del.AIC))]
factors = c('T', 'U', 'V', 'S', 'G', 'N')
T.valid = grep('T', names(mods2))
U.valid = grep('U', names(mods2))
V.valid = grep('V', names(mods2))
S.valid = grep('S', names(mods2))
G.valid = grep('G', names(mods2))
N.valid = grep('N', names(mods2))
plot.new()
plot.window(xlim = c(0.5, 7.5), ylim = c(0.5, length(del.AIC)+.5))
axis(2, at = 1:length(mods2), labels = gsub(pattern = 'lk.(\\w+)',replacement = "\\1", names(mods2)), las = 2)
axis(3, at = 1:7, labels = c(factors, expression(paste(Delta, 'AIC', sep = ''))))
for(ii in T.valid){ rr(c(1,ii))}
for(ii in U.valid){ rr(c(2,ii))}
for(ii in V.valid){ rr(c(3,ii))}
for(ii in S.valid){ rr(c(4,ii))}
for(ii in G.valid){ rr(c(5,ii))}
for(ii in N.valid){ rr(c(6,ii))}
for(ii in 1:length(del.AIC)){
  text(7, ii, paste(round(rev(del.AIC)[ii], 2)))
}
dev.off()



