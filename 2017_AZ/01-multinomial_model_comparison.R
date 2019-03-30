## Import data and run multinomila model comparison
## Fit to AZ data, from all seasons
## Fit a single age curve to all data
## Test three impriting hypothese: HA group level, HA subtype level, NA subtype level
## SAVE MODEL FITS 

## Clear memory
rm(list = ls())
#setwd('2017_AZ/')


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


#######################################
## SAVE MODEL FITS AND AIC
######################################
save(del.AIC, lk.A, lk.AG, lk.AN, lk.AS, file = modelfits)




## Re-do fits using bet-fit age pars for 2009 pandemic data

## Group-level imprinting
pandemic_AG = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AG$par[-c(1,2)], wPro.H1 = prog1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009, method = 'L-BFGS-B', lower = 0.001, upper = 1)

pandemic_AN = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AN$par[-c(1,2)], wPro.H1 = proN1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009, method = 'L-BFGS-B', lower = 0.001, upper = 1)

pandemic_AS = optim(par = c(rPro.H1 = .9), fn = nll_pandemic, fitted.age.pars = lk.AS$par[-c(1,2)], wPro.H1 = proH1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009, method = 'L-BFGS-B', lower = 0.001, upper = 1)

pandemic_A = list(value = nll_pandemic(par = c(rPro.H1 = 1), fitted.age.pars = lk.A$par, wPro.H1 = 1, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

## Get the difference between imprinting protection pars estimated for seasonal data and pandemic data
lk.AG$par['rPro.H1']-pandemic_AG$par
lk.AS$par['rPro.H1']-pandemic_AS$par
lk.AN$par['rPro.H1']-pandemic_AN$par

## Get the difference in likelihood between models fitting using pandemic and seasonal protection par
pandemic_AG$value
seasonal_AG = list(value = nll_pandemic(pars = lk.AG$par['rPro.H1'], fitted.age.pars = lk.AG$par, wPro.H1 = prog1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

pandemic_AN$value
seasonal_AN = list(value = nll_pandemic(pars = lk.AN$par['rPro.H1'], fitted.age.pars = lk.AN$par, wPro.H1 = proN1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

pandemic_AS$value
seasonal_AS = list(value = nll_pandemic(pars = lk.AS$par['rPro.H1'], fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009, dat.H1 = H1.master_2009, a0.4 = a0.4_2009, a5.10 = a5.10_2009, a11.17 = a11.17_2009, a18.24 = a18.24_2009, a25.31 = a25.31_2009, a32.38 = a32.38_2009, a39.45 = a39.45_2009, a46.52 = a46.52_2009, a53.59 = a53.59_2009, a60.66 = a60.66_2009, a67.73 = a67.73_2009, a74.80 = a74.80_2009, a81.90plus = a81.90plus_2009))

## Compare likelihoods 
seasonal_AG$value-pandemic_AG$value
seasonal_AS$value-pandemic_AS$value
seasonal_AN$value-pandemic_AN$value

## Compare estimates
AZ_seasonal_ests = c(AG = lk.AG$par['rPro.H1'], AN = lk.AN$par['rPro.H1'], AS = lk.AS$par['rPro.H1'])
AZ_pandemic_ests = c(AG = pandemic_AG$par['rPro.H1'], AN = pandemic_AN$par['rPro.H1'], AS = pandemic_AS$par['rPro.H1'])

## Compare fits
pdm_mods = mget(ls(pattern = "pandemic_A"))
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

save(AZ_seasonal_ests, AZ_pandemic_ests, pdm.del.AIC, pandemic_A, pandemic_AG, pandemic_AN, pandemic_AS, file = pdmfits)



### Comparge age distributions from seasons with different levels of drift
