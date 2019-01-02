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





