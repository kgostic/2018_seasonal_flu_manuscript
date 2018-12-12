## Calculate likelihood profiles for all paramters in models fitted to INSIGHT data
rm(list = ls())
source('00-Import_FLU002_-for-multinomial.R')
source('0func-multinomial_likelihood.R')
load('processed-data/fitted_models.RData')
attach(fits)


#############################
## Write a function wrapper to calculate likelihood profiles
############################
pvec = lk.A$par
lows = rep(.001, length(pvec))
highs = rep(5, length(pvec))

## Write a function wrapper to calculate one profile point value
## Input pvec, a named vector of best paramter estimates from a fitted likelihood. Should include the paramter to be fixed
##       fixpar, name of the parameter to be fixed
##       fixed.par.vale - value of paramter to be fixed
##       lows -  vector of lower par limits, corresponding to each entry in pvec, including the fixed parameter
##       highs - vector of upper par limits
one_prof_point = function(pvec, fixpar, fixed.par.value, lows, highs){
  drop = which(names(pvec) == fixpar) # Figure out index of fixed par
  # Optimize the likelihood with respect to all free paramters, other than the fixed par
  optim(par = pvec[-drop], fn = profile_function, fixed.par.name = fixpar, fixed.par.value = fixed.par.value, wPro.H1 = 1, dat.H1 = H1.master, wPro.H3 = 1, dat.H3 = H3.master, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, method = 'L-BFGS-B', lower = lows[-drop], upper = highs[-drop])$value
}