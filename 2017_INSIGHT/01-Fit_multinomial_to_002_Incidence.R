## Model comparison for INSIGHT FLU002
##  Fit multinomial models containing all possible combinations of imprinting hypotheses and medical history
##  Use AIC to compare model fits

#######################################
## Setup and import data
######################################
rm(list = ls())
source('00-Import_FLU002_-for-multinomial.R') ## Model inputs
source('0func-multinomial_likelihood.R') ## Likelihood function


## OUTPUTS
outfile1 = 'processed-data/fitted_models.RData'

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



pdf('INSIGHT_AIC.pdf', width = 4)
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
pdf('INSIGHT_AIC2.pdf', width = 4)
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
plotmod = function(pars, pro.H1 = 1, pro.H3 = 1, i.type = NULL){
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
  r74.80 = pars['r74.80'] 
  r81.90 = pars['r81.90'] 
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.risk = (b*r18.24*a18.24+b*r25.31*a25.31+ b*a32.38+ b*r39.45*a39.45+ b*r46.52*a46.52+ b*r53.59*a53.59+ b*r60.66*a60.66+ b*r67.73*a67.73+ b*r74.80*a74.80+ b*r81.90*a81.90)
  age.risk = age.risk/rowSums(age.risk)
  
  antivirals.H1 = (av.master*rAV+(1-av.master))
  antivirals.H3 = (av.master*rAV+(1-av.master))
  antivirals = rbind(antivirals.H1, antivirals.H3)
  
  underlying.H1 = (dx.master*rDX+(1-dx.master))
  underlying.H3 = (dx.master*rDX+(1-dx.master))
  underlying = rbind(underlying.H1, underlying.H3)
  
  vaccination.H1 = (vac.master*rVX.H1+(1-vac.master))
  vaccination.H3 = (vac.master*rVX.H3+(1-vac.master))
  
  imprinting.H1 = (pro.H1*rPro.H1+(1-pro.H1))
  imprinting.H3 = (pro.H3*rPro.H3+(1-pro.H3))
  
  par(mfrow = c(3,2))
  ## All rows in age baseline are the same, so plot one arbitrarily
  plot(18:90, age.risk[1,], main = 'Age effects', ylab = 'proportion of cases', xlab = 'case age')
  
  ## Antivial effects
  if(is.na(pars['rAV'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Antiviral treatment')
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(antivirals), col = 'black', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Antiviral treatment')
  }
  
  ## Underlying symptoms effects
  if(is.na(pars['rDX'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Underlying symptoms')
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(underlying), col = 'black', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Underlying symptoms')
  }
  
  ## Vaccination
  if(is.na(pars['rVX.H1'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Vaccination')
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(vaccination.H1), col = 'dodgerblue', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Vaccination')
    points(18:90, colMeans(vaccination.H3), col = 'firebrick1', cex = .7)
  }
  
  
  ## Imprinting protection
  if(is.na(pars['rPro.H1'])){
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste('Imprinting protection\n', i.type, sep = ''))
    text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(imprinting.H1), col = 'dodgerblue', cex = .7, ylab = 'est. relative risk impact', xlab = 'case age', main = 'Imprinting', ylim = c(0,1))
    points(18:90, colMeans(imprinting.H3), col = 'firebrick1', cex = .7)
  }
  
  
  ## Relative risk
  parnames = c("ages 18-24",  "ages 25-31", "fixed baseline",  "ages 39-45",  "ages 46-52",  "ages 53-59",  "ages 60-66", "ages 67-73",  "ages 74-80",  "ages 81-90", "", "", "", "antivirals", "underlying", "vaccination", paste('imprinting', i.type, sep = ', ')) 
  xvals = c(pars[grep(pattern = "r\\d\\d.\\d\\d", x = names(pars))], NA, NA, NA, pars['rAV'],pars['rDX'], pars['rVX.H1'], pars['rVX.H3'], pars['rPro.H1'], pars['rPro.H3'], 1)
  yvals = c(1, 2, 4:13, 14, 15, 16, 16, 17, 17, 3)
  par(mar = c(4, 7, 2, 1))
  plot(xvals, yvals, xlim = c(0, 1.3), xaxt = 'n', yaxt = 'n', xlab = 'Relative risk estimate', ylab = '', col = c(rep('black', 14), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), pch = 9)
  axis(side = 1, at = seq(0, 1.25, by = .25))
  axis(side = 2, at = 1:17, labels = parnames, las = 2)
  abline(v = 1, lty = 2)
  
  
  
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = tested.master/rowSums(tested.master)*age.risk * antivirals.H1 * underlying.H1 * vaccination.H1 * imprinting.H1
  pp.H3 = tested.master/rowSums(tested.master)* age.risk * antivirals.H3 * underlying.H3 * vaccination.H3 * imprinting.H3
  
  return(rbind(colSums(pp.H1/rowSums(pp.H1)*rowSums(H1.master)), colSums(pp.H3/rowSums(pp.H3)*rowSums(H3.master))))
}


## Plot best models
TS = plotmod(lk.TS$par, pro.H1 = proH1.master, pro.H3 = proH3.master, i.type = 'HA Subtype')

TN = plotmod(lk.TN$par, pro.H1 = proN1.master, pro.H3 = proN2.master, i.type = 'NA Subtype')

TT = plotmod(lk.T$par, pro.H1 = proH1.master, pro.H3 = proH3.master, i.type = 'HA Subtype')


pdf('INSIGHT002_predictions.pdf', width = 7, height = 4.5)
par(mfrow = c(1,2))
plot(18:90, colSums(H1.master), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H1N1')
lines(18:90, TS[1,], col = 'magenta', lwd = 2)
lines(18:90, TN[1,], col = 'orange', lwd = 2)
lines(18:90, TT[1,], col = 'dodgerblue', lwd = 2)
legend('topright', legend = c('TN', 'TS', 'T'), col = c('orange', 'magenta', 'dodgerblue'), lty = 1, bty = 'n')

plot(18:90, colSums(H3.master), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H3N2')
lines(18:90, TS[2,], col = 'magenta', lwd = 2)
lines(18:90, TN[2,], col = 'orange', lwd = 2)
lines(18:90, TT[2,], col = 'dodgerblue', lwd = 2)
dev.off()






