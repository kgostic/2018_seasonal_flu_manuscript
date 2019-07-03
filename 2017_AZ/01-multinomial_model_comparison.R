## Import data and run multinomila model comparison
## Fit to AZ data, from all seasons
## Fit a single age curve to all data
## Test three impriting hypothese: HA group level, HA subtype level, NA subtype level
## SAVE MODEL FITS 

## Clear memory
rm(list = ls())
setwd('~/Dropbox/R/2018_seasonal_flu/2017_AZ/')


## OUTPUTS
modelfits = 'processed-data/AZ_model_fits.RData'
pdmfits = 'processed-data/AZ_pandemic_fits.RData'
# outfile1 = '../figures/AZ_predictions.pdf' # Plot of model predictions vs. observed data
# outfile2 = '../figures/AZ_AIC.pdf' # Plot of AIC scores and factors included
# outfile3 = '../figures/AZ_NA_model_results.pdf' # Plot of AIC scores and factors included
# outfile4 = '../figures/AZ_HAsub_model_results.pdf' # Plot of AIC scores and factors included

#######################################
## Load data, model inputs, and likelihood function
######################################
source('00-Inputs_multinomial.R')
source('0func-likelihood.R')



## Test likelihood optimization
## Maximal model, includes imprinting protection
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
## 1. G 
lk.AG = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.AG

## 2. S
lk.AS = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.AS

## 3. N
lk.AN = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.AN

## 4. NULL
lk.A = nll.wrapper(pars.in = NULL, pro.H1 = 0, pro.H3 = 0, lower.in = NULL, upper = NULL); lk.A


# 
# ## 5. Age-specific to H1N1
# ## Concatenate vector to initialize pars that may or may not be included, with a vector of all age pars (always included)
# pvec = c(r5.10 = 1.1, r11.17 = .9, r18.24 = .9, r25.31 = .9, r32.38 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73= .9, r74.80 = .9, r81.90p = .9)
# A.H1 = optim(par = pvec, fn = nll_H1only, wPro.H1 = 0, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, dat.H1 = H1.master, method = 'L-BFGS-B', lower = c(rep(.001, 12)), upper = c(rep(5, 12))); A.H1
# 
# A.H3 = optim(par = pvec, fn = nll_H3only, wPro.H3 = 0, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, dat.H3 = H3.master, method = 'L-BFGS-B', lower = c(rep(.001, 12)), upper = c(rep(5, 12))); A.H3
# 
# 
# A.H1$value+A.H3$value



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


## Calculate akaike weights
raw.weights = exp(-.5*del.AIC)
AIC.weights = raw.weights/sum(raw.weights)
AIC.weights

## Summarize weights by imprinting type
imp.type = sub(pattern = "lk.\\w+?([NGS?])", replacement = "\\1", names(AIC.weights))
imp.type[grep('lk.', imp.type)] = 'none'

## Figure out how much weight each type gets
imp.type.weights = sapply(c('N', 'S', 'G', 'none'), FUN = function(tt){
  valid = AIC.weights[imp.type == tt]
  sum(valid)
})


AZ.imp.weights = imp.type.weights


#######################################
## SAVE MODEL FITS AND AIC
######################################
save(del.AIC, lk.A, lk.AG, lk.AN, lk.AS, AZ.imp.weights, file = modelfits)




## Re-do fits using bet-fit age pars for 2009 pandemic data

## Group-level imprinting
valid.row = 1 ## First, look at the 2008-09 season, roughly the first wave

firstwave_AG = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AG$par[-c(1,2)], wPro.H1 = prog1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,], method = 'L-BFGS-B', lower = 0.001, upper = 1)

firstwave_AN = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AN$par[-c(1,2)], wPro.H1 = proN1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,], method = 'L-BFGS-B', lower = 0.001, upper = 1)

firstwave_AS = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AS$par[-c(1,2)], wPro.H1 = proH1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,], method = 'L-BFGS-B', lower = 0.001, upper = 1)

firstwave_A = list(value = nll_pandemic(par = c(rPro.H1 = 1), fitted.age.pars = lk.AG$par[-c(1,2)], wPro.H1 = proH1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,]))

## Compare fits
pdm_mods = mget(ls(pattern = "firstwave_A"))
nll = numeric(length(pdm_mods))
pdmAICs = numeric(length(pdm_mods))
for(ii in 1:length(pdm_mods)){
  nll[ii] = pdm_mods[[ii]]$value
  pdmAICs[ii] = 2*length(pdm_mods[[ii]]$par)+2*pdm_mods[[ii]]$value
}

names(pdmAICs) = names(pdm_mods)
names(nll) = names(pdm_mods)
pdmAICs = sort(pdmAICs)
pdm.del.AIC = pdmAICs - min(pdmAICs)
pdm.del.AIC


## Get the difference between imprinting protection pars estimated for seasonal data and pandemic data
lk.AG$par['rPro.H1']-firstwave_AG$par
lk.AS$par['rPro.H1']-firstwave_AS$par
lk.AN$par['rPro.H1']-firstwave_AN$par

## Get the difference in likelihood between models fitting using pandemic and seasonal protection par
firstwave_AG$value
seasonal_AG = list(value = nll_pandemic(pars = lk.AG$par['rPro.H1'], fitted.age.pars = lk.AG$par, wPro.H1 = prog1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

firstwave_AN$value
seasonal_AN = list(value = nll_pandemic(pars = lk.AN$par['rPro.H1'], fitted.age.pars = lk.AN$par, wPro.H1 = proN1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

firstwave_AS$value
seasonal_AS = list(value = nll_pandemic(pars = lk.AS$par['rPro.H1'], fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

## Compare likelihoods 
seasonal_AG$value-firstwave_AG$value
seasonal_AS$value-firstwave_AS$value
seasonal_AN$value-firstwave_AN$value

## Compare estimates
AZ_seasonal_ests = c(AG = lk.AG$par['rPro.H1'], AN = lk.AN$par['rPro.H1'], AS = lk.AS$par['rPro.H1'])
AZ_firstwave_ests = c(AG = firstwave_AG$par['rPro.H1'], AN = firstwave_AN$par['rPro.H1'], AS = firstwave_AS$par['rPro.H1'])

## Compare fits
firstwavemods = mget(ls(pattern = "firstwave_A"))
nll = numeric(length(firstwavemods))
firstwaveAICs = numeric(length(firstwavemods))
for(ii in 1:length(firstwavemods)){
  nll[ii] = firstwavemods[[ii]]$value
  firstwaveAICs[ii] = 2*length(firstwavemods[[ii]]$par)+2*firstwavemods[[ii]]$value
}

names(firstwaveAICs) = names(firstwavemods)
names(nll) = names(firstwavemods)
firstwaveAICs = sort(firstwaveAICs)
firstwave.del.AIC = firstwaveAICs - min(firstwaveAICs)
firstwave.del.AIC






## repeate for second wave
## Group-level imprinting
valid.row = 2 

secondwave_AG = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AG$par[-c(1,2)], wPro.H1 = prog1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,], method = 'L-BFGS-B', lower = 0.001, upper = 1)

secondwave_AN = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AN$par[-c(1,2)], wPro.H1 = proN1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,], method = 'L-BFGS-B', lower = 0.001, upper = 1)

secondwave_AS = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AS$par[-c(1,2)], wPro.H1 = proH1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,], method = 'L-BFGS-B', lower = 0.001, upper = 1)

secondwave_A = list(value = nll_pandemic(par = c(rPro.H1 = 1), fitted.age.pars = lk.AG$par[-c(1,2)], wPro.H1 = proH1.master_2009[valid.row ,], dat.H1 = H1.master_2009[valid.row ,], a0.4 = a0.4_2009[valid.row,], a5.10 = a5.10_2009[valid.row,], a11.17 = a11.17_2009[valid.row,], a18.24 = a18.24_2009[valid.row,], a25.31 = a25.31_2009[valid.row,], a32.38 = a32.38_2009[valid.row,], a39.45 = a39.45_2009[valid.row,], a46.52 = a46.52_2009[valid.row,], a53.59 = a53.59_2009[valid.row,], a60.66 = a60.66_2009[valid.row,], a67.73 = a67.73_2009[valid.row,], a74.80 = a74.80_2009[valid.row,], a81.90plus = a81.90plus_2009[valid.row,]))

## Compare fits
pdm_mods = mget(ls(pattern = "secondwave_A"))
nll = numeric(length(pdm_mods))
pdmAICs = numeric(length(pdm_mods))
for(ii in 1:length(pdm_mods)){
  nll[ii] = pdm_mods[[ii]]$value
  pdmAICs[ii] = 2*length(pdm_mods[[ii]]$par)+2*pdm_mods[[ii]]$value
}

names(pdmAICs) = names(pdm_mods)
names(nll) = names(pdm_mods)
pdmAICs = sort(pdmAICs)
pdm.del.AIC = pdmAICs - min(pdmAICs)
pdm.del.AIC


## Get the difference between imprinting protection pars estimated for seasonal data and pandemic data
lk.AG$par['rPro.H1']-secondwave_AG$par
lk.AS$par['rPro.H1']-secondwave_AS$par
lk.AN$par['rPro.H1']-secondwave_AN$par

## Get the difference in likelihood between models fitting using pandemic and seasonal protection par
secondwave_AG$value
seasonal_AG = list(value = nll_pandemic(pars = lk.AG$par['rPro.H1'], fitted.age.pars = lk.AG$par, wPro.H1 = prog1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

secondwave_AN$value
seasonal_AN = list(value = nll_pandemic(pars = lk.AN$par['rPro.H1'], fitted.age.pars = lk.AN$par, wPro.H1 = proN1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

secondwave_AS$value
seasonal_AS = list(value = nll_pandemic(pars = lk.AS$par['rPro.H1'], fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

## Compare likelihoods 
seasonal_AG$value-secondwave_AG$value
seasonal_AS$value-secondwave_AS$value
seasonal_AN$value-secondwave_AN$value

## Compare estimates

AZ_secondwave_ests = c(AG = secondwave_AG$par['rPro.H1'], AN = secondwave_AN$par['rPro.H1'], AS = secondwave_AS$par['rPro.H1'])

## Compare fits
secondwavemods = mget(ls(pattern = "secondwave_A"))
nll = numeric(length(secondwavemods))
secondwaveAICs = numeric(length(secondwavemods))
for(ii in 1:length(secondwavemods)){
  nll[ii] = secondwavemods[[ii]]$value
  secondwaveAICs[ii] = 2*length(secondwavemods[[ii]]$par)+2*secondwavemods[[ii]]$value
}

names(secondwaveAICs) = names(secondwavemods)
names(nll) = names(secondwavemods)
secondwaveAICs = sort(secondwaveAICs)
secondwave.del.AIC = secondwaveAICs - min(secondwaveAICs)
secondwave.del.AIC

save(AZ_seasonal_ests, AZ_firstwave_ests, firstwave.del.AIC, firstwave_A, firstwave_AG, firstwave_AN, firstwave_AS, AZ_secondwave_ests, secondwave.del.AIC, secondwave_A, secondwave_AG, secondwave_AN, secondwave_AS, file = pdmfits)



### Comparge age distributions from seasons with different levels of drift
