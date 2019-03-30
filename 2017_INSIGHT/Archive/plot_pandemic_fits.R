### Plot pandemic model fitss
setwd('../../2018_seasonal_flu/')
rm(list = ls())

## Load INSIGHT pandemic and seasonal fits
load('2017_INSIGHT/processed-data/INSIGHT_fitted_models.RData')
load('2017_INSIGHT/processed-data/INSIGHT_fitted_pandemic_models.RData')

# Load AZ pandemic and seasonal fits
load('2017_AZ/processed-data/AZ_model_fits.RData')
load('2017_AZ/processed-data/AZ_pandemic_fits.RData')

## Load AZ inputs
setwd('2017_AZ/')
source('00-Inputs_multinomial.R')
setwd('../')

## Set up panel layout
layout(matrix(c(1,2,3,3), nrow = 2, ncol = 2))


## Define a function with which to calculate AZ pandemic fits
plotmodAZ = function(rPro.H1, seasonal_pars, pro.H1){
  
  b = 1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r5.10 = seasonal_pars['r5.10'] # Expected risk for 5 to 10 year olds
  r11.17 = seasonal_pars['r11.17'] # etc.
  r18.24 = seasonal_pars['r18.24'] 
  r25.31 = seasonal_pars['r25.31']
  r32.38 = seasonal_pars['r32.38']
  r39.45 = seasonal_pars['r39.45']
  r46.52 = seasonal_pars['r46.52'] 
  r53.59 = seasonal_pars['r53.59']
  r60.66 = seasonal_pars['r60.66'] 
  r67.73 = seasonal_pars['r67.73'] 
  r74.80 = seasonal_pars['r74.80'] 
  r81.90p = seasonal_pars['r81.90p'] 
  
  
  ## Calculate predicted case distributions, as in the likelihood ##
  age.baseline = (b*a0.4_2009 +b*r5.10*a5.10_2009 +b*r11.17*a11.17_2009+b*r18.24*a18.24_2009+b*r25.31*a25.31_2009+ b*r32.38*a32.38_2009+ b*r39.45*a39.45_2009+ b*r46.52*a46.52_2009+ b*r53.59*a53.59_2009+ b*r60.66*a60.66_2009+ b*r67.73*a67.73_2009+ b*r74.80*a74.80_2009+ b*r81.90p*a81.90plus_2009)
  age.baseline.rr = age.baseline # Save unnormalized version
  age.baseline = age.baseline/sum(age.baseline)
  
  
  
  imprinting.H1 = (pro.H1*rPro.H1+(1-pro.H1))
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = age.baseline * imprinting.H1
  ## Return the predicted age distributions of infection. Plot these below against the observed data
  return((pp.H1/sum(pp.H1)*sum(H1.master_2009))[as.character(1918:2009)])
}
# Get prediction for model AS, but don't save to pdf
spred = plotmodAZ(rPro.H1 = lk.AS$par['rPro.H1'], seasonal_pars = lk.AS$par, pro.H1 = proH1.master_2009)
spred_pdm =  plotmodAZ(rPro.H1 = pandemic_AS$par, seasonal_pars = lk.AS$par, pro.H1 = proH1.master_2009)
npred = plotmodAZ(rPro.H1 = lk.AN$par['rPro.H1'], seasonal_pars = lk.AN$par, pro.H1 = proN1.master_2009)
npred_pdm =  plotmodAZ(rPro.H1 = pandemic_AN$par, seasonal_pars = lk.AN$par, pro.H1 = proN1.master_2009)


## Plot panel A, AZ 2009 pandemic data and model fits
plot(-1918:-2009, H1.master_2009[as.character(1918:2009)], type = 'l', lwd = 2, xlab = 'birth year', ylab = 'case count', main = 'Model fits: H1N1', xaxt = 'n')
axis(side = 1, at = seq(-2009, -1918, by = 10), labels = seq(2009, 1918, by = -10), line = 0)
lines(-1918:-2009, spred, col = 'green', lwd = 1)
lines(-1918:-2009, spred_pdm, col = 'blue', lwd = 1)
plot(-1918:-2009, H1.master_2009[as.character(1918:2009)], type = 'l', lwd = 2, xlab = 'birth year', ylab = 'case count', main = 'Model fits: H1N1', xaxt = 'n')
axis(side = 1, at = seq(-2009, -1918, by = 10), labels = seq(2009, 1918, by = -10), line = 0)
lines(-1918:-2009, npred, col = 'green', lwd = 2)
lines(-1918:-2009, npred_pdm, col = 'blue', lwd = 2)
legend('topright', legend = c('observed', 'AN fit', 'AS fit', 'AG fit', 'A   fit'), col = c('black', cols[5:2]), lwd = c(2,2,2,1,1), bty = 'n')
mtext('E', side = 3, line = 1, at = -2020, font = 2)