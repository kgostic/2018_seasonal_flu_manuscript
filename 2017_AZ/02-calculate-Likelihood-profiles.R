## Calculate likelihood profiles for all model fits
########################################
# SETUP
########################################
## Clear memory
rm(list = ls())
## Load model fits and source model inputs and likelihood functions
load('processed-data/AZ_model_fits.RData') # Model fits saved as lk.A, lk.AG, lk.AS, lk.AN
source('0func-likelihood.R')
source('00-Inputs_multinomial.R')
## Remove unnecessary vars
rm('bys', 'ii', 'validH1', 'validH3', 'weights')
library('parallel')


## OUTPUTS
outfile1 = 'processed-data/AZ_risk_profiles.RData'



########################################
# Write a function wrapper to calculate profiles
########################################
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
optim(par = pvec[-drop], fn = profile_func, fixed.par.name = fixpar, fixed.par.value = fixed.par.value, wPro.H1 = 1, dat.H1 = H1.master, wPro.H3 = 1, dat.H3 = H3.master, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, method = 'L-BFGS-B', lower = lows[-drop], upper = highs[-drop])$value
}

# ## Test
# one_prof_point(pvec = lk.AG$par, fixpar = 'r5.10', fixed.par.value = .5, lows = rep(.001, length(lk.AG$par)), highs = c(1,1, rep(5, 12)))
# 
# one_prof_point(pvec = lk.AG$par, fixpar = 'rPro.H1', fixed.par.value = .5, lows = rep(.001, length(lk.AG$par)), highs = c(1,1, rep(5, 12)))



###########################
## Profiles for age-specific risk
###########################
###### Initialize storage
grid = seq(.005, 2, by = .005) # Define grid of relative risk points to test for each paramter
# Store provile neg log likelihood values in a matrix, with the fixed parameter listed on rows, and grid value listed on columns
A.age.prof = AG.age.prof = AS.age.prof = AN.age.prof = matrix(NA, nrow = 12, ncol = length(grid), dimnames = list(names(lk.A$par), grid))

## Run a profile grid for each age paramter
## Run grids in parallel
cl = makeCluster(detectCores()-1) # Make cluster
clusterExport(cl, ls()) # Export all variables to cluster
age.pars = names(lk.A$par)

## For each age paramter
for(pp in age.pars){
## lk.A age profiles
A.age.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.A$par, fixpar = pp, lows = rep(.001, 12), highs = rep(5, 12))

AG.age.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.AG$par, fixpar = pp, lows = rep(.001, 14), highs = c(1,1, rep(5, 12)))

AS.age.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.AS$par, fixpar = pp, lows = rep(.001, 14), highs = c(1,1, rep(5, 12)))

AN.age.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.AN$par, fixpar = pp, lows = rep(.001, 14), highs = c(1,1, rep(5, 12)))
}
stopCluster(cl)





###########################
## Profiles for imprinting protection pars
###########################
###### Initialize storage
grid = seq(.005, 1, by = .005) # Define grid of relative risk points to test for each paramter
# Store provile neg log likelihood values in a matrix, with the fixed parameter listed on rows, and grid value listed on columns
AG.imp.prof = AS.imp.prof = AN.imp.prof = matrix(NA, nrow = 2, ncol = length(grid), dimnames = list(c('rPro.H1', 'rPro.H3'), grid))
cl = makeCluster(detectCores()-1) # Make cluster
clusterExport(cl, ls()) # Export all variables to cluster
pro.pars = c('rPro.H1', 'rPro.H3')
## For each age paramter
for(pp in pro.pars){
  AG.imp.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.AG$par, fixpar = pp, lows = rep(.001, 14), highs = c(1,1, rep(5, 12)))
  
  AS.imp.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.AS$par, fixpar = pp, lows = rep(.001, 14), highs = c(1,1, rep(5, 12)))
  
  AN.imp.prof[pp, ] = parSapply(cl = cl, X = grid, FUN = one_prof_point, pvec = lk.AN$par, fixpar = pp, lows = rep(.001, 14), highs = c(1,1, rep(5, 12)))
}
stopCluster(cl)

save(A.age.prof, AS.age.prof, AN.age.prof, AG.age.prof, AS.imp.prof, AG.imp.prof, AN.imp.prof, file = outfile1)




## Get profile CIs

# Define a function to calculate the likelihood ratio threshold, as a function of the best nll value, and the number of constrained pars (df = 1 if only one par is fixed in the profile)
LR.Threshold = function(NLL_best, df){
  #-2log(LR) is distributed chi2 with df given by the number of fixed parameters in the profile
  #algebraically, the threshold for being in the CI is given by:
  threshold = qchisq(.95, df)/2+NLL_best
  threshold
}

# Define a function to extract the CI bounds from the grid, and the profile values
LR.CI = function(threshold, nll.vec, pars.vec){
  if(any(which(nll.vec > threshold) < which.min(nll.vec))  &  any(which(nll.vec > threshold) > which.min(nll.vec)) ){ #If the minimum is not an end point
    #Find string before and after the min value
    lower = nll.vec[1:which.min(nll.vec)]
    upper = nll.vec[(which.min(nll.vec)-1):length(nll.vec)]
    #Extract low CI from first string, and upper CI from upper string
    CI = c(pars.vec[which.min(abs(lower-threshold))], pars.vec[length(lower)+which.min(abs(upper-threshold))])
  }else{
    #If the first value is the minimum
    if(any(which(nll.vec > threshold) < which.min(nll.vec)) == FALSE){ CI = c(pars.vec[1], pars.vec[which.min(abs(nll.vec-threshold))]) }
    #If the last value is the maximum
    if(any(which(nll.vec > threshold) > which.min(nll.vec)) == FALSE){ CI = c(pars.vec[which.min(abs(nll.vec-threshold))], pars.vec[length(nll.vec)])}
  }
  CI
}





      