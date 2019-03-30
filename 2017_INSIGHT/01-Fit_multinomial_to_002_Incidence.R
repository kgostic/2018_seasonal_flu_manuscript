## Model comparison for INSIGHT FLU002
##  Fit multinomial models containing all possible combinations of imprinting hypotheses and medical history
##  Use AIC to compare model fits

#######################################
## Setup and import data
######################################
rm(list = ls())
setwd('~/Dropbox/R/2018_seasonal_flu/2017_INSIGHT/')
source('00-Import_FLU002_-for-multinomial.R') ## Model inputs
source('0func-multinomial_likelihood.R') ## Likelihood function


## OUTPUTS
outfile1 = 'processed-data/INSIGHT_fitted_models.RData'
outfile2 = 'processed-data/INSIGHT_fitted_pandemic_models.RData'
outfile3 = 'processed-data/INSIGHT_pandemic_and_seasonal_imp_risk_ests.RData'

## Test likelihood optimization
# # Maximal model
## pars.in gives a vector of initial guesses for free paramters. Must be a nmed vector.
## pro.H1 is the master matrix describing probabilities of protection against H1N1. proH1.master describes HA subtye-specific protection, prog1.master describes HA group-specific protection, or proN1.master describes N1-specific protection
## pro.H3 - master protection matrix relevant to H3N2 protection
## lower.in - vector of lower limits for estimated free par values. Order: rAv, rDX, rVS, rpro.H1, rpro.H3
## upper.in - upper limits for estimated free par values
nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(.001, .001, 1, 1, .001, .001, .001, .001), upper = c(3, 3, 10, 10, 1, 1, 1, 1))


## You can exclude factors from the model by excluding their paramter names from the pars.in vector. This signals the wrapper to substitute the null relative risk value (1), which effectively removes the factor from the model.
## Reduced model, vaccination and imprinting only
nll.wrapper(pars.in = c('rVX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(.001, .001, .001, .001), upper = c(1, 1, 1, 1))





#######################################
## Model comparison 
##      null is all cases tested
##      fit to H1N1
######################################
## set upper and lower bounds for all parameter estiamtes
av.low = .001; av.high = 10 # Antiviral treatment could reasonably be associated with increased risk (given that confirmed cases are more likely to be treated), or decreased risk. Allow values above and below 1.
vx.low = .001; vx.high = 1 ## Vaccination should predict decreased risk. Prohibit valeus above 1.
dx.low = .001; dx.high = 10 # Underlying symptoms -- similar to AV treatment
pro.low = 0.001; pro.high = 1 ## Protection should only decrease risk. Prohibit values above 1.

## All models contain age, so A abbreviation omitted.
# T = antivirals
# U = underlying symptoms
# V = vaccination
# G = group
# S = subtype
# N = neuraminidase

## 1. TUVG 
lk.ATUVG = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, dx.high, vx.high, vx.high, pro.high, pro.high)); lk.ATUVG

## 2. TUVS
lk.ATUVS = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, dx.high, vx.high, vx.high, pro.high, pro.high)); lk.ATUVS

## 3. TUVN
lk.ATUVN = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, dx.high, vx.high, vx.high, pro.high, pro.high)); lk.ATUVN

## 4. TUV
lk.ATUV = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(av.low, dx.low, vx.low, vx.low), upper = c(av.high, dx.high, vx.high, vx.high)); lk.ATUV

## 5. TUG 
lk.ATUG = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, dx.low, pro.low, pro.low), upper = c(av.high, dx.high, pro.high, pro.high)); lk.ATUG

## 6. TUS
lk.ATUS = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, dx.low, pro.low, pro.low), upper = c(av.high, dx.high, pro.high, pro.high)); lk.ATUS

## 7. TUN
lk.ATUN = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, dx.low, pro.low, pro.low), upper = c(av.high, dx.high, pro.high, pro.high)); lk.ATUN

## 8. TU
lk.ATU = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(av.low, dx.low), upper = c(av.high, dx.high)); lk.ATU

## 9. TVG 
lk.ATVG = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, vx.high, vx.high, pro.high, pro.high)); lk.ATVG

## 10. TVS
lk.ATVS = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, vx.high, vx.high, pro.high, pro.high)); lk.ATVS

## 11. TVN
lk.ATVN = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, vx.high, vx.high, pro.high, pro.high)); lk.ATVN

## 12. TV
lk.ATV = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(av.low, vx.low, vx.low), upper = c(av.high, vx.high, vx.high)); lk.ATV

## 13. UVG 
lk.AUVG = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(dx.high, vx.high, vx.high, pro.high, pro.high)); lk.AUVG

## 14. UVS
lk.AUVS = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(dx.high, vx.high, vx.high, pro.high, pro.high)); lk.AUVS

## 15. UVN
lk.AUVN = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(dx.high, vx.high, vx.high, pro.high, pro.high)); lk.AUVN

## 16. UV
lk.AUV = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(dx.low, vx.low, vx.low), upper = c(dx.high, vx.high, vx.high)); lk.AUV

## 17. UG 
lk.AUG = nll.wrapper(pars.in = c('rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(dx.low, pro.low, pro.low), upper = c(dx.high, pro.high, pro.high)); lk.AUG

## 18. US
lk.AUS = nll.wrapper(pars.in = c('rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(dx.low, pro.low, pro.low), upper = c(dx.high, pro.high, pro.high)); lk.AUS

## 19. UN
lk.AUN = nll.wrapper(pars.in = c('rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(dx.low, pro.low, pro.low), upper = c(dx.high, pro.high, pro.high)); lk.AUN

## 20. U
lk.AU = nll.wrapper(pars.in = c('rDX' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(dx.low), upper = c(dx.high)); lk.AU

## 21. TG 
lk.ATG = nll.wrapper(pars.in = c('rAV' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, pro.low, pro.low), upper = c(av.high, pro.high, pro.high)); lk.ATG

## 22. TS
lk.ATS = nll.wrapper(pars.in = c('rAV' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, pro.low, pro.low), upper = c(av.high, pro.high, pro.high)); lk.ATS

## 23. TN
lk.ATN = nll.wrapper(pars.in = c('rAV' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, pro.low, pro.low), upper = c(av.high, pro.high, pro.high)); lk.ATN

## 24. T
lk.AT = nll.wrapper(pars.in = c('rAV' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(av.low), upper = c(av.high)); lk.AT

## 25. VG 
lk.AVG = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(vx.low, vx.low, pro.low, pro.low), upper = c(vx.high, vx.high, pro.high, pro.high)); lk.AVG

## 26. VS
lk.AVS = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(vx.low, vx.low, pro.low, pro.low), upper = c(vx.high, vx.high, pro.high, pro.high)); lk.AVS

## 27. VN
lk.AVN = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(vx.low, vx.low, pro.low, pro.low), upper = c(vx.high, vx.high, pro.high, pro.high)); lk.AVN

## 28. V
lk.AV = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = 1, pro.H3 = 1, lower.in = c(vx.low, vx.low), upper = c(vx.high, vx.high)); lk.AV

## 29. G 
lk.AG = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.AG

## 30. S
lk.AS = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.AS

## 31. N
lk.AN = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.AN

## 32. NULL
lk.A = nll.wrapper(pars.in = NULL, pro.H1 = 1, pro.H3 = 1, lower.in = NULL, upper = NULL); lk.A






## Pull out the variable names that store likelihoods
mods = mget(ls(pattern = "lk."))
nlls = numeric(length(mods))
AICs = numeric(length(mods))
for(ii in 1:length(mods)){
  nlls[ii] = mods[[ii]]$value
  AICs[ii] = 2*length(mods[[ii]]$par)+2*mods[[ii]]$value
}

names(AICs) = names(mods)
names(nlls) = names(mods)
AICs = sort(AICs)
del.AIC = AICs - min(AICs)
del.AIC.full = del.AIC
del.AIC
# T = antivirals
# U = underlying symptoms
# V = vaccination
# G = group
# S = subtype
# N = neuraminidase



#####################
## save fitted models
#####################
fits = mget(ls(pattern = 'lk.\\w+'))
save(fits, del.AIC, file = outfile1)










#####################
## re-fit imprinting protection during 2009 pandemic, using best-fit age curves from above
#####################
rm(pdm.indices)
## 1. TUVG 
pdm.ATUVG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATUVG$par, pro.H1 = prog1.master_pandemic); pdm.ATUVG

## 2. TUVS
pdm.ATUVS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATUVS$par, pro.H1 = proH1.master_pandemic); pdm.ATUVS

## 3. TUVN
pdm.ATUVN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATUVN$par, pro.H1 = proN1.master_pandemic); pdm.ATUVN

## 4. TUV
pdm.ATUV = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.ATUV$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.ATUV

## 4. TUG 
pdm.ATUG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATUG$par, pro.H1 = prog1.master_pandemic); pdm.ATUG

## 5. TUS
pdm.ATUS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATUS$par, pro.H1 = proH1.master_pandemic); pdm.ATUS

## 6. TUN
pdm.ATUN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATUN$par, pro.H1 = proN1.master_pandemic); pdm.ATUN

## 7. TU
pdm.ATU = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.ATU$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.ATU


## 9. TVG 
pdm.ATVG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATVG$par, pro.H1 = prog1.master_pandemic); pdm.ATVG

## 10. TVS
pdm.ATVS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATVS$par, pro.H1 = proH1.master_pandemic); pdm.ATVS

## 11. TVN
pdm.ATVN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATVN$par, pro.H1 = proN1.master_pandemic); pdm.ATVN

## 12. TV
pdm.ATV = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.ATV$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = 1, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.ATV

## 13. UVG 
pdm.AUVG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AUVG$par, pro.H1 = prog1.master_pandemic); pdm.AUVG

## 14. UVS
pdm.AUVS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AUVS$par, pro.H1 = proH1.master_pandemic); pdm.AUVS

## 15. UVN
pdm.AUVN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AUVN$par, pro.H1 = proN1.master_pandemic); pdm.AUVN

## 16. UV
pdm.AUV = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.AUV$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.AUV

## 17. UG 
pdm.AUG = nll.wrapper.pdm(pars.in = c( 'rPro.H1' = .5), seasonal_fit = lk.AUG$par, pro.H1 = prog1.master_pandemic); pdm.AUG

## 18. US
pdm.AUS = nll.wrapper.pdm(pars.in = c( 'rPro.H1' = .5), seasonal_fit = lk.AUS$par, pro.H1 = proH1.master_pandemic); pdm.AUS

## 19. UN
pdm.AUN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AUN$par, pro.H1 = proN1.master_pandemic); pdm.AUN

## 20. U
pdm.AU = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.AU$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.AU

## 21. TG 
pdm.ATG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATG$par, pro.H1 = prog1.master_pandemic); pdm.ATG

## 22. TS
pdm.ATS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATS$par, pro.H1 = proH1.master_pandemic); pdm.ATS

## 23. TN
pdm.ATN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.ATN$par, pro.H1 = proN1.master_pandemic); pdm.ATN

## 24. T
pdm.AT = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.AT$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.AT

## 25. VG 
pdm.AVG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AVG$par, pro.H1 = prog1.master_pandemic); pdm.AVG

## 26. VS
pdm.AVS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AVS$par, pro.H1 = proH1.master_pandemic); pdm.AVS

## 27. VN
pdm.AVN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AVN$par, pro.H1 = proN1.master_pandemic); pdm.AVN

## 28. V
pdm.AV = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.AV$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.AV

## 29. G 
pdm.AG = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AG$par, pro.H1 = prog1.master_pandemic); pdm.AG

## 30. S
pdm.AS = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AS$par, pro.H1 = proH1.master_pandemic); pdm.AS

## 31. N
pdm.AN = nll.wrapper.pdm(pars.in = c('rPro.H1' = .5), seasonal_fit = lk.AN$par, pro.H1 = proN1.master_pandemic); pdm.AN

## 32. A
pdm.A = list(value = nll_pandemic(pars = c(rPro.H1 =1), fitted.age.pars = lk.A$par, wAV = av.master_pandemic, wDX = dx.master_pandemic, wVX = vac.master_pandemic, wPro.H1 = prog1.master_pandemic, dat.H1 = H1.master_pandemic, a18.24 = a18.24_pandemic, a25.31 = a25.31_pandemic, a32.38 = a32.38_pandemic, a39.45 = a39.45_pandemic, a46.52 = a46.52_pandemic, a53.59 = a53.59_pandemic, a60.66 = a60.66_pandemic, a67.73 = a67.73_pandemic, a74.80 = a74.80_pandemic, a81.90 = a81.90_pandemic, tested.master = tested.master_pandemic/rowSums(tested.master_pandemic))); pdm.A






## Pull out the variable names that store likelihoods
pdmmods = mget(ls(pattern = "pdm.A"))
nlls = numeric(length(pdmmods))
pdmAICs = numeric(length(pdmmods))
for(ii in 1:length(pdmmods)){
  nlls[ii] = pdmmods[[ii]]$value
  pdmAICs[ii] = 2*length(pdmmods[[ii]]$par)+2*pdmmods[[ii]]$value
}

names(pdmAICs) = names(pdmmods)
names(nlls) = names(pdmmods)
pdmAICs = sort(pdmAICs)
del.pdmAIC = pdmAICs - min(pdmAICs)
del.pdmAIC.full = del.pdmAIC
del.pdmAIC




#####################
## save fitted models
#####################
pdm_mods = pdmmods
save(pdm_mods, del.pdmAIC, file = outfile2)





#####################
## compare imprinting protection levels from models fitted to seasonal data vs. pandemic data
#####################
seasonal_ests = sapply(fits, function(xx) xx$par['rPro.H1'])
pandemic_ests = sapply(pdm_mods, function(xx) xx$par['rPro.H1'], simplify = TRUE)
pandemic_ests[which( sapply(pandemic_ests, length) == 0)] = NA

seasonal_ests = seasonal_ests[-which(is.na(seasonal_ests))]
pandemic_ests = pandemic_ests[-which(is.na(pandemic_ests))]
pandemic_ests = unlist(pandemic_ests)

## higher protection coefficient means weaker protection. Assume pandemic will be higher and subtract seasonal - pandemic. If this yields a positive value, then pandemic protection is stronger
seasonal_ests-pandemic_ests

plot((seasonal_ests-pandemic_ests))
gl = grep(names(seasonal_ests), pattern = 'G.')
points(gl, (seasonal_ests[gl]-pandemic_ests[gl]), col = 'red')
sl = grep(names(seasonal_ests), pattern = 'S.')
points(sl, (seasonal_ests-pandemic_ests)[sl], col = 'blue')
nl = grep(names(seasonal_ests), pattern = 'N.')
points(nl, (seasonal_ests-pandemic_ests)[nl], col = 'green')
abline(h = 0)


save(seasonal_ests, pandemic_ests, file = outfile3)
