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
## Note: a lower bound of exactly 0 will occasionally return a NAN value and crash the optimizer.
##       use a positive, lower bound close to 0 (e.g. 0.001) to avoid this
##      choose a lower bound closer to 0 if you MLEs are crashing to 0.001
## 1. G 
lk.AG = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), # Inital guess
                    pro.H1 = prog1.master, pro.H3 = prog2.master, # Input fixed imprinting probabilities
                    lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); # Set lower bounds on parameter estimates
lk.AG # Return

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


