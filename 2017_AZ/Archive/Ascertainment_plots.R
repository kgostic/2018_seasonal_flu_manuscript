rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')
load('H1ests_vv_OAS_2017-09-03.RData')
load('H3ests_vv_OAS2017-09-03.RData')


## Try to account for ascertainment
# Data
cols = gray.colors(10)[-1]
plot(0:99, H3.inputs[1, ]/sum(H3.inputs[1, ]), ylim = c(0, .09), main = 'Imprinting', col = cols[1])
for(ii in 1:nrow(H3.inputs)){ points(0:99, H3.inputs[ii, ]/sum(H3.inputs[ii, ]), col = cols[ii])}
points(0:99, colSums(H3.inputs)/sum(H3.inputs), pch = 16)

###### ------ No Imprinting ------ ######
# Simulate model best fit using OAS estimates
OAS = simulate.incidence(alpha.hat = H3.ests.vv.OAS$par['alpha.hat'], AA.hat = H3.ests.vv.OAS$par['AA.hat'], R0.hat = H3.ests.vv.OAS$par['R0.hat'], vv.hat = H3.ests.vv.OAS$par['vv.hat'], tau.hat = H3.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)
# Plot using a blue line
use.years = rownames(H3.inputs)
lines(0:99, colSums(OAS[use.years, ])/sum(OAS[use.years, ]), col = 'black', lwd = 2)


###### ------ Ascertainment: Age ------ ######
# Simulate model best fit using best estimates
load('H3ests_vv_OAS_age2017-09-04.RData')
raw.age = simulate.incidence(alpha.hat = H3.ests.vv.OAS.age$par['alpha.hat'], AA.hat = H3.ests.vv.OAS.age$par['AA.hat'], R0.hat = H3.ests.vv.OAS.age$par['R0.hat'], vv.hat = H3.ests.vv.OAS.age$par['vv.hat'], tau.hat = H3.ests.vv.OAS.age$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)

# Modulate the raw simulation based on case ascertainment parameters
elderly.mat = H3.inputs*0; elderly.mat[,as.character(65:99)] = 1
adult.child.mat = 1-elderly.mat
age =  raw.age[use.years, ] * (elderly.mat * H3.ests.vv.OAS.age$par['age.par'] + adult.child.mat)
lines(0:99, colSums(age)/sum(age), col = 'red4', lwd = 2, lty = 2)


###### ------ Ascertainment: Imprinting  ------ ######
# Simulate model best fit using best estimates
load('H3ests_vv_OAS_imp2017-09-03.RData')
raw.imp = simulate.incidence(alpha.hat = H3.ests.vv.OAS.imp$par['alpha.hat'], AA.hat = H3.ests.vv.OAS.imp$par['AA.hat'], R0.hat = H3.ests.vv.OAS.imp$par['R0.hat'], vv.hat = H3.ests.vv.OAS.imp$par['vv.hat'], tau.hat = H3.ests.vv.OAS.imp$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)

# Modulate the raw simulation based on case ascertainment parameters
imp =  raw.imp[use.years, ] * (wm.H3[use.years, ] + H3.ests.vv.OAS.imp$par['HH']*wo.H3[use.years, ])
lines(0:99, colSums(imp)/sum(imp), col = 'blue3', lwd = 2, lty = 2)




###### ------ Ascertainment: Imprinting x Age ------ ######
# Simulate model best fit using best estimates
load('H3ests_vv_OAS_age_imp2017-09-03.RData')
raw.age.imp = simulate.incidence(alpha.hat = H3.ests.vv.OAS.age.imp$par['alpha.hat'], AA.hat = H3.ests.vv.OAS.age.imp$par['AA.hat'], R0.hat = H3.ests.vv.OAS.age.imp$par['R0.hat'], vv.hat = H3.ests.vv.OAS.age.imp$par['vv.hat'], tau.hat = H3.ests.vv.OAS.age.imp$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)
# Modulate the raw simulation based on case ascertainment parameters
elderly.mat = H3.inputs*0; elderly.mat[,as.character(65:99)] = 1
adult.child.mat = 1-elderly.mat
age.imp =  raw.age.imp[use.years, ] * (wm.H3[use.years, ] + H3.ests.vv.OAS.age.imp$par['HH']*wo.H3[use.years, ]) * (elderly.mat * H3.ests.vv.OAS.age.imp$par['age.par'] + adult.child.mat)
lines(0:99, colSums(age.imp)/sum(age.imp), col = 'purple4', lwd = 2)


###### ------ Ascertainment:  Age sloped ------ ######
# Simulate model best fit using best estimates
load('H3ests_vv_OAS_age_slope2017-09-04.RData')
raw.age.slope = simulate.incidence(alpha.hat =  H3.ests.vv.OAS.age_slope$par['alpha.hat'], AA.hat =  H3.ests.vv.OAS.age_slope$par['AA.hat'], R0.hat =  H3.ests.vv.OAS.age_slope$par['R0.hat'], vv.hat =  H3.ests.vv.OAS.age_slope$par['vv.hat'], tau.hat =  H3.ests.vv.OAS.age_slope$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)
# Modulate the raw simulation based on case ascertainment parameters
elderly.mat = matrix(c(rep(0, 65), 1:35*H3.ests.vv.OAS.age_slope$par['age.par']/35+1), nrow = nrow(H3.inputs), ncol = 100, byrow = T); colnames(elderly.mat) = 0:99
adult.child.mat = elderly.mat*0; adult.child.mat[, as.character(0:64)] = 1
age.slope =  raw.age.slope[use.years, ] * (elderly.mat + adult.child.mat)
lines(0:99, colSums(age.slope)/sum(age.slope), col = 'red', lwd = 2)



###### ------ Ascertainment: Imprinting x Age sloped ------ ######
# Simulate model best fit using best estimates
load('H3ests_vv_OAS_age_slope_imp2017-09-03.RData')
raw.age.imp = simulate.incidence(alpha.hat = H3.ests.vv.OAS.age_slope.imp$par['alpha.hat'], AA.hat = H3.ests.vv.OAS.age_slope.imp$par['AA.hat'], R0.hat = H3.ests.vv.OAS.age_slope.imp$par['R0.hat'], vv.hat = H3.ests.vv.OAS.age_slope.imp$par['vv.hat'], tau.hat = H3.ests.vv.OAS.age_slope.imp$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)
# Modulate the raw simulation based on case ascertainment parameters
elderly.mat = matrix(c(rep(0, 65), 1:35*20/35+1), nrow = nrow(H3.inputs), ncol = 100, byrow = T); colnames(elderly.mat) = 0:99
adult.child.mat = elderly.mat*0; adult.child.mat[, as.character(0:64)] = 1
age.imp =  raw.age.imp[use.years, ] * (wm.H3[use.years, ] + 9*wo.H3[use.years, ]) * (elderly.mat + adult.child.mat)
lines(0:99, colSums(age.imp)/sum(age.imp), col = 'blue', lwd = 2)



elderly.mat = matrix(c(rep(0, 65), 1:35*10/35+1), nrow = nrow(H3.inputs), ncol = 100, byrow = T); colnames(elderly.mat) = 0:99
adult.child.mat = elderly.mat*0; adult.child.mat[, as.character(0:64)] = 1
xx = (elderly.mat * 3 + adult.child.mat)

# ## Child age imprinting
# load('H3ests_vv_OAS_age_imp_child2017-09-04.RData')
# raw.age.imp.child = simulate.incidence(alpha.hat = H3.ests.vv.OAS.age.imp.child$par['alpha.hat'], AA.hat = H3.ests.vv.OAS.age.imp.child$par['AA.hat'], R0.hat = H3.ests.vv.OAS.age.imp.child$par['R0.hat'], vv.hat = H3.ests.vv.OAS.age.imp.child$par['vv.hat'], tau.hat = H3.ests.vv.OAS.age.imp.child$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)

# elderly.mat = H3.inputs*0; elderly.mat[,as.character(65:99)] = 1
# adult.child.mat = 1-elderly.mat
# child.mat = c(rep(1, 5), rep(0,95))
# adult.child.mat[1:5] = 0
# age.imp.child = raw.age.imp.child[use.years, ] * (elderly.mat *  H3.ests.vv.OAS.age.imp.child$par['age.par'] + adult.child.mat + H3.ests.vv.OAS.age.imp.child$par['child.par']*child.mat )
# lines(0:99, colSums(age.imp.child)/sum(age.imp.child), col = 'green')


legend('topright', 
       c('Baseline', 'Age', 'Imprinting', 'Age x imprinting', 'Age sloped', 'Age sloped x imprinting'), 
       lty = c(1, 2, 2, 1, 2, 1), col = c('black', 'red4', 'blue3', 'purple4', 'red', 'purple'), bty = 'n')
