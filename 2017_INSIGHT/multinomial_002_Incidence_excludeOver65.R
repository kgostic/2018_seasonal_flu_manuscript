## Model comparison for INSIGHT FLU002


#######################################
## Setup and import data
######################################
rm(list = ls())
setwd('~/Dropbox/R/2017_INSIGHT/')
source('Import_FLU002_multinomial.R')



#######################################
## Define likelihood
######################################
#### INPUTS:
####    pars - named vector of pars to be fit. If not all four pars exist in a given model, code automatically assigns the null value of 1. Parameters not named in "pars" will not be fit.
####    null - vector, or matrix describing null age distribution. May be given by all cases tested, or all confirmed flu cases.
####    wAV - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_i) who used antivirals
####    wDX - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_i) who had underlying symptoms
####    wVX - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_i) who were vaccinated in the previous 12 months
####    wPro - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_i) who are protected against the challenge strain based on childhood imprinting
####    xx - n-vector or mxn matrix of observed case counts. Each entry represents the number of cases of a given age (n_i), observed in a given country and season (m_i).

##### OUTPUTS:
#####    negative log likelihood (scalar)


##### NOTES:  
####    All models tested are nested in this likelihood structure. E.g. to exclude the effects of protection from the model, input wPro = 0, and omit "rPro" from pars.

nll = function(pars, wAV.H1, wDX.H1, wVX.H1, wPro.H1, dat.H1, wAV.H3, wDX.H3, wVX.H3, wPro.H3, tested.master, dat.H3, a18.24, a25.31, a32.38, a39.45, a46.52, a53.59, a60.66, a67.73){
  # 1. Assign parameters to be fit
  rAV = ifelse(is.na(pars['rAV']), 1, pars['rAV']) # Relative risk given antiviral use
  rDX = ifelse(is.na(pars['rDX']), 1, pars['rDX']) # Relative risk given underlying symptoms
  rVX.H1 = ifelse(is.na(pars['rVX.H1']), 1, pars['rVX.H1']) # Relative risk given vaccination
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rVX.H3 = ifelse(is.na(pars['rVX.H3']), 1, pars['rVX.H3']) # Relative risk given vaccination
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  r18.24 = pars['r18.24'] # Baseline expectation for 18-24 year olds
  r25.31 = pars['r25.31']
  b =1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73']
 
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.risk = (b*r18.24*a18.24+b*r25.31*a25.31+ b*a32.38+ b*r39.45*a39.45+ b*r46.52*a46.52+ b*r53.59*a53.59+ b*r60.66*a60.66+ b*r67.73*a67.73)
  age.risk = age.risk/rowSums(age.risk)

 
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = tested.master*age.risk * (wAV.H1*rAV+(1-wAV.H1))* (wDX.H1*rDX+(1-wDX.H1))* (wVX.H1*rVX.H1+(1-wVX.H1))* (wPro.H1*rPro.H1+(1-wPro.H1))
  
  pp.H3 = tested.master*age.risk * (wAV.H3*rAV+(1-wAV.H3))* (wDX.H3*rDX+(1-wDX.H3))* (wVX.H3*rVX.H3+(1-wVX.H3))* (wPro.H3*rPro.H3+(1-wPro.H3))
  
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
use = as.character(18:73)

# # Test function
nll(pars = c(rAV = .5, rDX = .5, rVX.H1 = .5, rPro.H1 = .5, rVX.H3 = .5, rPro.H3 = .5, r18.24 = .9, r25.31 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73 = .9), tested.master = (tested.master/rowSums(tested.master))[,use], wAV.H1 = av.master[,use], wAV.H3 = av.master[,use], wDX.H1 = dx.master[,use], wDX.H3 = dx.master[,use], wVX.H1 = vac.master[,use], wVX.H3 = vac.master[,use], wPro.H1 = proH1.master[,use], wPro.H3 = proH3.master[,use], a18.24 = a18.24[,use], a25.31 = a25.31[,use], a32.38 = a32.38[,use], a39.45 = a39.45[,use], a46.52 = a46.52[,use], a53.59 = a53.59[,use], a60.66 = a60.66[,use], a67.73 = a67.73[,use], dat.H1 = H1.master[,use], dat.H3 = H3.master[,use])


## Write a function wrapper, so that you only have to input a vector of initial par values and par value limits involved in model comparison
nll.wrapper = function(pars.in, pro.H1, pro.H3, lower.in, upper.in){
  ## Concatenate vector to initialize pars that may or may not be included, with a vector of all age pars (always included)
  pvec = c(pars.in, r18.24 = .9, r25.31 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73 = .9)
  ex = as.character(65:90)
 optim(par = pvec, fn = nll, wAV.H1 = av.master[,use], wAV.H3 = av.master[,use], wDX.H1 = dx.master[,use], wDX.H3 = dx.master[,use], wVX.H1 = vac.master[,use], wVX.H3 = vac.master[,use], wPro.H1 = pro.H1[,use], wPro.H3 = pro.H3[,use], tested.master = (tested.master[,use]/rowSums(tested.master[,use])), a18.24 = a18.24[,use], a25.31 = a25.31[,use], a32.38 = a32.38[,use], a39.45 = a39.45[,use], a46.52 = a46.52[,use], a53.59 = a53.59[,use], a60.66 = a60.66[,use], a67.73 = a67.73[,use], dat.H1 = H1.master[,use], dat.H3 = H3.master[,use], method = 'L-BFGS-B', lower = c(lower.in, rep(.001, 7)), upper = c(upper.in, rep(5, 7)))
}

# Test
# Maximal model
nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(.001, .001, 1, 1, .001, .001, .001, .001), upper = c(3, 3, 10, 10, 1, 1, 1, 1))

## Reduced model, vaccination and imprinting only
nll.wrapper(pars.in = c('rVX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(.001, .001, .001, .001), upper = c(1, 1, 1, 1))



#######################################
## Model comparison 
##      null is all cases tested
##      fit to H1N1
######################################
## set upper and lower bounds for all pars
av.low = .001; av.high = 10
vx.low = .001; vx.high = 1
dx.low = 1; dx.high = 10
pro.low = 0.001; pro.high = 1

filler_1s = (H1.master*0)[,use]+1

# T = antivirals
# U = underlying symptoms
# V = vaccination
# G = group
# S = subtype
# N = neuraminidase
## 1. TUVG 
lk.TUVG = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, dx.high, vx.high, vx.high, pro.high, pro.high)); lk.TUVG

## 2. TUVS
lk.TUVS = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, dx.high, vx.high, vx.high, pro.high, pro.high)); lk.TUVS

## 3. TUVN
lk.TUVN = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, dx.high, vx.high, vx.high, pro.high, pro.high)); lk.TUVN

## 4. TUV
lk.TUV = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(av.low, dx.low, vx.low, vx.low), upper = c(av.high, dx.high, vx.high, vx.high)); lk.TUV

## 5. TUG 
lk.TUG = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, dx.low, pro.low, pro.low), upper = c(av.high, dx.high, pro.high, pro.high)); lk.TUG

## 6. TUS
lk.TUS = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, dx.low, pro.low, pro.low), upper = c(av.high, dx.high, pro.high, pro.high)); lk.TUS

## 7. TUN
lk.TUN = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, dx.low, pro.low, pro.low), upper = c(av.high, dx.high, pro.high, pro.high)); lk.TUN

## 8. TU
lk.TU = nll.wrapper(pars.in = c('rAV' = .5, 'rDX' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(av.low, dx.low, vx.low, vx.low), upper = c(av.high, dx.high)); lk.TU

## 9. TVG 
lk.TVG = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, vx.high, vx.high, pro.high, pro.high)); lk.TVG

## 10. TVS
lk.TVS = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, vx.high, vx.high, pro.high, pro.high)); lk.TVS

## 11. TVN
lk.TVN = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, vx.low, vx.low, pro.low, pro.low), upper = c(av.high, vx.high, vx.high, pro.high, pro.high)); lk.TVN

## 12. TV
lk.TV = nll.wrapper(pars.in = c('rAV' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(av.low, vx.low, vx.low), upper = c(av.high, vx.high, vx.high)); lk.TV

## 13. UTG 
lk.UTG = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(dx.high, vx.high, vx.high, pro.high, pro.high)); lk.UTG

## 14. UTS
lk.UTS = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(dx.high, vx.high, vx.high, pro.high, pro.high)); lk.UTS

## 15. UTN
lk.UTN = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(dx.low, vx.low, vx.low, pro.low, pro.low), upper = c(dx.high, vx.high, vx.high, pro.high, pro.high)); lk.UTN

## 16. UT
lk.UT = nll.wrapper(pars.in = c('rDX' = .5, 'rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(dx.low, vx.low, vx.low), upper = c(dx.high, vx.high, vx.high)); lk.UT

## 17. UG 
lk.UG = nll.wrapper(pars.in = c('rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(dx.low, pro.low, pro.low), upper = c(dx.high, pro.high, pro.high)); lk.UG

## 18. US
lk.US = nll.wrapper(pars.in = c('rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(dx.low, pro.low, pro.low), upper = c(dx.high, pro.high, pro.high)); lk.US

## 19. UN
lk.UN = nll.wrapper(pars.in = c('rDX' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(dx.low, pro.low, pro.low), upper = c(dx.high, pro.high, pro.high)); lk.UN

## 20. U
lk.U = nll.wrapper(pars.in = c('rDX' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(dx.low), upper = c(dx.high)); lk.U

## 21. TG 
lk.TG = nll.wrapper(pars.in = c('rAV' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(av.low, pro.low, pro.low), upper = c(av.high, pro.high, pro.high)); lk.TG

## 22. TS
lk.TS = nll.wrapper(pars.in = c('rAV' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(av.low, pro.low, pro.low), upper = c(av.high, pro.high, pro.high)); lk.TS

## 23. TN
lk.TN = nll.wrapper(pars.in = c('rAV' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(av.low, pro.low, pro.low), upper = c(av.high, pro.high, pro.high)); lk.TN

## 24. T
lk.T = nll.wrapper(pars.in = c('rAV' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(av.low), upper = c(av.high)); lk.T

## 25. VG 
lk.VG = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(vx.low, vx.low, pro.low, pro.low), upper = c(vx.high, vx.high, pro.high, pro.high)); lk.VG

## 26. VS
lk.VS = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(vx.low, vx.low, pro.low, pro.low), upper = c(vx.high, vx.high, pro.high, pro.high)); lk.VS

## 27. VN
lk.VN = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5, 'rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(vx.low, vx.low, pro.low, pro.low), upper = c(vx.high, vx.high, pro.high, pro.high)); lk.VN

## 28. V
lk.V = nll.wrapper(pars.in = c('rVX.H1' = .5, 'rVX.H3' = .5), pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = c(vx.low, vx.low), upper = c(vx.high, vx.high)); lk.V

## 29. G 
lk.G = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = prog1.master, pro.H3 = prog2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.G

## 30. S
lk.S = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proH1.master, pro.H3 = proH3.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.S

## 31. N
lk.N = nll.wrapper(pars.in = c('rPro.H1' = .5, 'rPro.H3' = .5), pro.H1 = proN1.master, pro.H3 = proN2.master, lower.in = c(pro.low, pro.low), upper = c(pro.high, pro.high)); lk.N

## 32. NULL
lk.null = nll.wrapper(pars.in = NULL, pro.H1 = filler_1s, pro.H3 = filler_1s, lower.in = NULL, upper = NULL); lk.null






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
del.AIC.full = del.AIC
del.AIC
# T = antivirals
# U = underlying symptoms
# V = vaccination
# G = group
# S = subtype
# N = neuraminidase


## Remove degenerate models
ii = 1
while(ii < length(del.AIC)){
 remainder = del.AIC[(ii+1):length(del.AIC)] - del.AIC[ii]
 drop = which(round(remainder, 3) %% 2 == 0)
 del.AIC = del.AIC[!names(del.AIC) %in% names(drop)]
 ii = ii+1
}
del.AIC



## Plot AIC results
mods = mods[rev(names(del.AIC.full))]
factors = c('T', 'U', 'V', 'S', 'G', 'N')
T.valid = grep('T', names(mods))
U.valid = grep('U', names(mods))
V.valid = grep('V', names(mods))
S.valid = grep('S', names(mods))
G.valid = grep('G', names(mods))
N.valid = grep('N', names(mods))
rr = function(center, width = 1, height = 1, col.in = 'navy'){
  rect(xleft = center[1]-width/2, ybottom = center[2]-height/2, xright = center[1]+width/2, ytop = center[2]+height/2, col = col.in, border = 'black')
}
dev.off(); 
pdf('INSIGHT_AICexover73.pdf', width = 4)
plot.new()
plot.window(xlim = c(0.5, 7.5), ylim = c(0.5, 32.5))
axis(2, at = 1:32, labels = gsub(pattern = 'lk.(\\w+)',replacement = "\\1", names(mods)), las = 2)
axis(3, at = 1:7, labels = c(factors, expression(paste(Delta, 'AIC', sep = ''))))
for(ii in T.valid){ rr(c(1,ii, width = .8))}
for(ii in U.valid){ rr(c(2,ii))}
for(ii in V.valid){ rr(c(3,ii))}
for(ii in S.valid){ rr(c(4,ii))}
for(ii in G.valid){ rr(c(5,ii))}
for(ii in N.valid){ rr(c(6,ii), width = .8)}
for(ii in 1:32){
  text(7, ii, paste(round(rev(del.AIC.full)[ii], 2)))
}

## Gray out degenerate models
degenerate = which(!rev(del.AIC.full) %in% del.AIC)
for(ii in degenerate){
  rr(center = c(3.45, ii), width = 5.8, height = 1, col.in = 'gray')
}

dev.off()
pdf('INSIGHT_AIC2exover73.pdf', width = 4)
layout(c(1,1,2,2,2))
plot.new()
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




## Plot model results
plotmod = function(pars, pro.H1 = filler_1s, pro.H3 = filler_1s, i.type = NULL){
  # 1. Assign parameters to be fit
  rAV = ifelse(is.na(pars['rAV']), 1, pars['rAV']) # Relative risk given antiviral use
  rDX = ifelse(is.na(pars['rDX']), 1, pars['rDX']) # Relative risk given underlying symptoms
  rVX.H1 = ifelse(is.na(pars['rVX.H1']), 1, pars['rVX.H1']) # Relative risk given vaccination
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rVX.H3 = ifelse(is.na(pars['rVX.H3']), 1, pars['rVX.H3']) # Relative risk given vaccination
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  r18.24 = pars['r18.24'] # Baseline expectation for 18-24 year olds
  r25.31 = pars['r25.31']
  b =1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73']
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.risk = (b*r18.24*a18.24[,use]+b*r25.31*a25.31[,use]+ b*a32.38[,use]+ b*r39.45*a39.45[,use]+ b*r46.52*a46.52[,use]+ b*r53.59*a53.59[,use]+ b*r60.66*a60.66[,use]+b*r67.73*a67.73[,use])
  age.risk = age.risk/rowSums(age.risk)
  
  antivirals.H1 = (av.master[,use]*rAV+(1-av.master[,use]))
  antivirals.H3 = (av.master[,use]*rAV+(1-av.master[,use]))
  antivirals = rbind(antivirals.H1, antivirals.H3)
  
  underlying.H1 = (dx.master[,use]*rDX+(1-dx.master[,use]))
  underlying.H3 = (dx.master[,use]*rDX+(1-dx.master[,use]))
  underlying = rbind(underlying.H1, underlying.H3)
  
  vaccination.H1 = (vac.master[,use]*rVX.H1+(1-vac.master[,use]))
  vaccination.H3 = (vac.master[,use]*rVX.H3+(1-vac.master[,use]))
  
  imprinting.H1 = (pro.H1[,use]*rPro.H1+(1-pro.H1[,use]))
  imprinting.H3 = (pro.H3[,use]*rPro.H3+(1-pro.H3[,use]))
  
  par(mfrow = c(3,2))
  ## All rows in age baseline are the same, so plot one arbitrarily
  plot(18:73, age.risk[1,], main = 'Age effects', ylab = 'proportion of cases', xlab = 'case age')
  
  ## Antivial effects
  if(is.na(pars['rAV'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Antiviral treatment')
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:73, colMeans(antivirals), col = 'black', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Antiviral treatment')
  }
  
  ## Underlying symptoms effects
  if(is.na(pars['rDX'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Underlying symptoms')
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:73, colMeans(underlying), col = 'black', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Underlying symptoms')
  }
  
  ## Vaccination
  if(is.na(pars['rVX.H1'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Vaccination')
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:73, colMeans(vaccination.H1), col = 'dodgerblue', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Vaccination')
    points(18:73, colMeans(vaccination.H3), col = 'firebrick1', cex = .7)
  }
  
  
  ## Imprinting protection
  if(is.na(pars['rPro.H1'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste('Imprinting protection\n', i.type, sep = ''))
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:73, colMeans(imprinting.H1), col = 'dodgerblue', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Imprinting', ylim = c(0,1))
    points(18:73, colMeans(imprinting.H3), col = 'firebrick1', cex = .7)
  }

  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = age.risk * antivirals.H1 * underlying.H1 * vaccination.H1 * imprinting.H1
  pp.H3 = age.risk * antivirals.H3 * underlying.H3 * vaccination.H3 * imprinting.H3
  
  return(rbind(colSums(pp.H1/rowSums(pp.H1)*rowSums(H1.master)), colSums(pp.H3/rowSums(pp.H3)*rowSums(H3.master))))
}


## Plot best models
TS = plotmod(lk.TS$par, pro.H1 = proH1.master[,use], pro.H3 = proH3.master[,use], i.type = 'HA Group')

TN = plotmod(lk.TN$par, pro.H1 = proN1.master[,use], pro.H3 = proN2.master[,use], i.type = 'NA Subtype')

TG = plotmod(lk.TG$par, pro.H1 = filler_1s, pro.H3 = filler_1s, i.type = NULL)


pdf('INSIGHT002_predictions_excludeover73.pdf', width = 7, height = 4.5)
par(mfrow = c(1,2))
plot(18:73, colSums(H1.master[,use]), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H1N1')
lines(18:73, TS[1,], col = 'magenta', lwd = 2)
lines(18:73, TN[1,], col = 'orange', lwd = 2)
lines(18:73, TG[1,], col = 'dodgerblue', lwd = 2)
legend('topright', legend = c('TS', 'TN', 'TG'), col = c('magenta', 'orange', 'dodgerblue'), lty = 1, bty = 'n')

plot(18:73, colSums(H3.master[,use]), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H3N2')
lines(18:73, TS[2,], col = 'magenta', lwd = 2)
lines(18:73, TN[2,], col = 'orange', lwd = 2)
lines(18:73, TG[2,], col = 'dodgerblue', lwd = 2)
dev.off()






