######################### Plot model fits ################################
## Run from seasonal_flu home directory
## Start with results from AZ data
rm(list = ls())
## setwd('~/Dropbox/R/2018_seasonal_flu') ## Run from the seasonal_flu home folder
## Load libraries
library(viridis)
library(scales)
## Temporarily switch to the 2017_AZ subdirectory and load AZ model fits, CIs and model inputs
## We will load the INISHGT inputs below, when we start that plot.
setwd('2017_AZ/')
load('processed-data/AZ_model_fits.RData')
load('processed-data/AZ_CIs.RData')
source('00-Inputs_multinomial.R') 
setwd('../')# Switch back into the home directory


## OUTPUTS
outfile1 = 'figures/AZ_model_fits.pdf' ## Output one figure that plots AZ model fits to data, AIC comparison, and individuals factors from the best model
outfile2 = 'figures/INSIGHT_model_fits.pdf' ## Output a second figure that plots the same results for INSIGHT data


## Set up color palette
cols = rev(viridis_pal(alpha = .7, begin = 1, end = .1)(5))
show_col(cols)


## Define a function analagous to the likelihood above, which outputs model predictions
#### INPUTS
##     pars - named vector of paramter values (should be best estimates from fits above)
##     pro.H1 - probabilities of H1N1 protection
##     pro.H3 - probabilities of H3N2 protection
##     i.type - character, type of imprinting protection (can be 'HA subtype', 'HA group', or 'NA subtype'). Used in plot titles.
plotmod1 = function(pars, CIs, pro.H1 = 1, pro.H3 = 1, i.type = NULL){
  
  ## Parse paramter inputs ##
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  b = 1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r5.10 = pars['r5.10'] # Expected risk for 5 to 10 year olds
  r11.17 = pars['r11.17'] # etc.
  r18.24 = pars['r18.24'] 
  r25.31 = pars['r25.31']
  r32.38 = pars['r32.38']
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73'] 
  r74.80 = pars['r74.80'] 
  r81.90p = pars['r81.90p'] 
  
  
  ## Calculate predicted case distributions, as in the likelihood ##
  age.baseline = (b*a0.4 +b*r5.10*a5.10 +b*r11.17*a11.17+b*r18.24*a18.24+b*r25.31*a25.31+ b*r32.38*a32.38+ b*r39.45*a39.45+ b*r46.52*a46.52+ b*r53.59*a53.59+ b*r60.66*a60.66+ b*r67.73*a67.73+ b*r74.80*a74.80+ b*r81.90p*a81.90plus)
  age.baseline.rr = age.baseline # Save unnormalized version
  age.baseline = age.baseline/rowSums(age.baseline)
  imprinting.H1 = (pro.H1*rPro.H1+(1-pro.H1))
  imprinting.H3 = (pro.H3*rPro.H3+(1-pro.H3))
  
  ## Age-specific risk prediction
  ## Plot row 15, which represents the last observed season, and aligns 0-year-olds with the first column (2015 birth year)
  ## No need to take the mean, because age-specific predictions are identical across countries and years
  par(mar = c(3.5, 4, 3.5, 5))
  plot(0:97, age.baseline.rr[15,], main = 'Age effects', ylab = 'relative risk', xlab = 'age', bty = 'n', ylim = c(0,1))
  mtext('A', side = 3, line = 1.5, at = -5, font = 2)
  abline(h = 1, lty = 2)
  
  ## Predicted effects from imprinting protection
  ## Plot the colmeans (mean predicted protection across birth years)
  ## This is necessary because imprinting reconstructions differ slightly from year to year, as children (with recent birth years) get older.
  if(is.na(pars['rPro.H1'])){ # If no imprinting protection par, this factor was not relevant to the model fit, and we plot an empty windwo
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste('Imprinting protection\n', i.type, sep = ''), bty = 'n')
    text(.5, .5, 'NA', cex = 2)
  }else{ # Else, if imprinting protection was included in the model, plot mean birth year-specific relative risk, after adjusting for imprinting protection.
    ## Plot just one row here, because predicted protection changes with birth year
    plot(-2015:-1918, imprinting.H1[15,], col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'birth year', main = 'Imprinting effects', ylim = c(0,1), bty = 'n', xaxt = 'n')
    axis(side = 1, at = seq(-2015, -1918, by = 10), labels = seq(2015, 1918, by = -10), line = 0)
    points(-2015:-1918, imprinting.H3[15,], col = 'firebrick1', cex = .7)
    }
  ## H1N1 risk in blue, H3N2 risk in red
  mtext('B', side = 3, line = 1.5, at = -5, font = 2)
  abline(h = 1, lty = 2)
  
  ## Plot relative risk point estimates for each free parameter.
  par(mar = par('mar')+c(.5,3,.5,.5))
  parnames = c("ages 0-4", "ages5-10", "ages 11-17", "ages 18-24",  "ages 25-31", "ages 32-38",  "ages 39-45",  "ages 46-52",  "ages 53-59",  "ages 60-66", "ages 67-73",  "ages 74-80",  "ages 81+", "", "", "", paste('Impr.', i.type, sep = ', ')) 
  xvals = c(NA, pars[grep(pattern = "r\\d\\d?.\\d\\d", x = names(pars))], NA, NA, NA, pars['rPro.H1'], pars['rPro.H3'])
  ## and CIS
  xlows = CIs[1,]
  xlows = c(NA, xlows[grep(pattern = "r\\d\\d?.\\d\\d", x = colnames(CIs))], NA, NA, NA, xlows['rPro.H1'], xlows['rPro.H3'])
  xhis = CIs[2,]
  xhis = c(NA, xhis[grep(pattern = "r\\d\\d?.\\d\\d", x = colnames(CIs))], NA, NA, NA, xhis['rPro.H1'], xhis['rPro.H3'])
  yvals = c(1:15, 16, 16.9, 17.1)
  #print(rbind(yvals, xvals, xlows, xhis, c(rep('black', 12), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black')))
  #par(mar = c(4, 7, 2, 1))
  plot(xvals, yvals, xlim = c(0, 1.6), xaxt = 'n', yaxt = 'n', xlab = 'relative risk', ylab = '', col = c(rep('black', 16), 'dodgerblue', 'firebrick1'), pch = 9, main = "Maximum likelihood estimates")
  segments(x0 = xlows, y0 = yvals, x1 = xhis, col = c(rep('black', 16), 'dodgerblue', 'firebrick1'), lwd = 3)
  axis(side = 1, at = seq(0, 1.6, by = .25))
  axis(side = 2, at = 1:17, labels = parnames, las = 2)
  abline(v = 1, lty = 2)
  mtext('C', side = 3, line = 1.5, at = -.35, font = 2, xpd = NA)
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = age.baseline * imprinting.H1
  pp.H3 = age.baseline * imprinting.H3
  ## Return the predicted age distributions of infection. Plot these below against the observed data
  return(rbind(rev(colSums(pp.H1/rowSums(pp.H1)*rowSums(H1.master))), rev(colSums(pp.H3/rowSums(pp.H3)*rowSums(H3.master)))))
}
# Get prediction for model AS, but don't save to pdf
spred = plotmod1(lk.AS$par, pro.H1 = proH1.master, pro.H3 = proH3.master, i.type = 'HA Sub', CIs = AS.CIs) # Get predicted age distributions using hte best fit subtype-specific imprinting model.
gpred = plotmod1(lk.AG$par, pro.H1 = prog1.master, pro.H3 = prog2.master, i.type = "HA grp", CIs = AG.CIs)
apred = plotmod1(lk.A$par, pro.H1 = 1, pro.H3 = 1, i.type = NULL, CIs = A.CIs)






{
###########################################################################
############################ START AZ PDF #################################
###########################################################################
## ------------------------------------------------------------------------
pdf(outfile1, height = 8)
layout(mat = matrix(c(1,3, 2,3, 4,4, 5,6), ncol = 2, byrow = T), widths = c(1,1), heights = c(1,1,1,1.4))
par(mgp = c(2,1,0))
#layout.show(6)
##########################
## 1-3
##########################
npred = plotmod1(lk.AN$par, pro.H1 = proN1.master, pro.H3 = proN2.master, i.type = 'NA Sub', CIs = AN.CIs) # Get predicted age distributions 

##########################
## 4. Plot AIC rankings
##########################
par(mar = par('mar')+c(0, 4, 0, 4))
## Write a function to plot a filled rectangle of fixed size, given the center coordinate
rr = function(center, width = 1, height = 1, col.in = 'navy'){
  rect(xleft = center[1]-width/2, ybottom = center[2]-height/2, xright = center[1]+width/2, ytop = center[2]+height/2, col = col.in, border = 'white')
}
## Repeat, only plotting non-degenerate models
## Plot AIC results
mods = mget(ls(pattern = 'lk.\\w+'))
mods2 = mods[rev(names(del.AIC))]
factors = c('A', 'S', 'G', 'N')
S.valid = grep('S', names(mods2))
G.valid = grep('G', names(mods2))
N.valid = grep('N', names(mods2))
A.valid = 1:length(mods2)
plot.new()
plot.window(xlim = c(0.5, 5.5), ylim = c(0.5, length(del.AIC)+.5))
mtext(text = 'Model comparison', side = 3, line = 2, font = 2)
segments(x1 = -.25, x0 = 5.5, y0 = 4.75, xpd = NA)
axis(2, at = 1:length(mods2), labels = gsub(pattern = 'lk.(\\w+)',replacement = "\\1", names(mods2)), las = 2)
axis(3, at = 1:5, labels = c(factors, expression(paste(Delta, 'AIC', sep = ''))), tick = FALSE, cex.axis = 1.2, font = 2, line = -.5)
mtext(text = 'Model', side = 2, line = 2.5, cex = .8)
for(ii in A.valid){ rr(c(1,ii), col.in = cols[ii+1])}
for(ii in S.valid){ rr(c(2,ii), col.in = cols[ii+1])}
for(ii in G.valid){ rr(c(3,ii), col.in = cols[ii+1])}
for(ii in N.valid){ rr(c(4,ii), col.in = cols[ii+1])}
for(ii in 1:length(del.AIC)){
  text(5, ii, paste(round(rev(del.AIC)[ii], 2)))
}
mtext('D', side = 3, line = 2, at = -.25, font = 2)


##########################
## 5-6. Plot model fits to observed data
##########################
par(mar = c(3,3,.5,2)+.5)
plot(-1918:-2015, rev(colSums(H1.master)), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'case count', main = 'Model fits: H1N1', xaxt = 'n')
axis(side = 1, at = seq(-2015, -1918, by = 10), labels = seq(2015, 1918, by = -10), line = 0)
lines(-1918:-2015, apred[1,], col = cols[2], lwd = 1)
lines(-1918:-2015, gpred[1,], col = cols[3], lwd = 1)
lines(-1918:-2015, spred[1,], col = cols[4], lwd = 2)
lines(-1918:-2015, npred[1,], col = cols[5], lwd = 2)
legend('topright', legend = c('observed', 'AN fit', 'AS fit', 'AG fit', 'A   fit'), col = c('black', cols[5:2]), lwd = c(2,2,2,1,1), bty = 'n')
mtext('E', side = 3, line = 1, at = -2020, font = 2)

plot(-1918:-2015, rev(colSums(H3.master)), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'case count', main = 'Model fits: H3N2', xaxt = 'n')
axis(side = 1, at = seq(-2015, -1918, by = 10), labels = seq(2015, 1918, by = -10), line = 0)
lines(-1918:-2015, apred[2,], col = cols[2], lwd = 1)
lines(-1918:-2015, gpred[2,], col = cols[3], lwd = 1)
lines(-1918:-2015, spred[2,], col = cols[4], lwd = 2)
lines(-1918:-2015, npred[2,], col = cols[5], lwd = 2)
mtext('F', side = 3, line = 1, at = -2020, font = 2)
dev.off()
## ------------------------------------------------------------------------
###########################################################################
###########################################################################
}









###########################################################################
## Plot INSIGHT results
###########################################################################
## Temporarily switch to the 2017_INSIGHT subdirectory and load INSIGHT model fits, CIs and model inputs
rm(list = c(ls(pattern = 'lk.\\w+')))
setwd('2017_INSIGHT/')
load('processed-data/INSIGHT_fitted_models.RData')
CIs = read.csv('processed-data/INSIGHT_CIs.csv')
colnames(CIs)[3:4] = c('low', 'high')
source('00-Import_FLU002_-for-multinomial.R') 
setwd('../')# Switch back into the home directory


# cols2 = viridis_pal(alpha = .7, begin = .95, end = 0, option = "A")(32-5)
# show_col(c(cols, cols2))
cols1 = rev(c(rev(cols),rep('gray', 32-5)))
show_col(cols1)

## Define a function analagous to the likelihood above, which outputs model predictions
#### INPUTS
##     pars - named vector of paramter values (should be best estimates from fits above)
##     pro.H1 - probabilities of H1N1 protection
##     pro.H3 - probabilities of H3N2 protection
##     i.type - character, type of imprinting protection (can be 'HA subtype', 'HA group', or 'NA subtype'). Used in plot titles.
plotmod = function(modname){
  mod = fits[[paste('lk.', modname, sep = "")]]
  pars = mod$par
  pro.H1 = 1
  pro.H3 = 1
  CIplot = CIs[CIs$modname == modname, ]
  i.type = NULL
  ## Set up protection inputs based on model name
  if(grepl('N', modname)){pro.H1 = proN1.master; pro.H3 = proN2.master; i.type = 'NA sub'} # If the model considers NA subtype-level protection, use NA-specific protection inputs
  if(grepl('S', modname)){pro.H1 = proH1.master; pro.H3 = proH3.master; i.type = 'HA sub'} # If the model considers HA subtype-level protection, use HA-specific protection inputs
  if(grepl('G', modname)){pro.H1 = prog1.master; pro.H3 = prog2.master; i.type = 'HA grp'} # If the model considers HA group-level protection, use group-specific protection inputs
  
  

  
  
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
  age.rr = age.risk # save non-normalized version for plot
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


  ## All rows in age baseline are the same, so plot one arbitrarily
  plot(18:90, age.rr[1,], main = 'Age effects', ylab = 'relative risk', xlab = 'age', ylim = c(.5, 1.5))
  abline(h = 1, lty = 2, bty = 'n')
  mtext(text = 'A', side = 3, line = 1, at = 10, font = 2)

  ## Antivial effects
  if(is.na(pars['rAV'])){
    # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Antiviral effects')
    # text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(antivirals), col = 'black', cex = .7, ylab = 'relative risk', xlab = 'age', main = 'Antiviral effects', ylim = c(.5, 1.5), bty = 'n')
    abline(h = 1, lty = 2)
  }
  mtext(text = 'B', side = 3, line = 1, at = 10, font = 2)

  ## Underlying symptoms effects
  if(is.na(pars['rDX'])){
    # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Underlying cond. effects')
    # text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(underlying), col = 'black', cex = .7, ylab = 'relative risk', xlab = 'age', main = 'Underlying cond. effects', ylim = c(.5, 1.5), bty = 'n')
    abline(h = 1, lty = 2)
    mtext(text = 'C', side = 3, line = 1, at = 10, font = 2)
  }
 

  ## Vaccination
  if(is.na(pars['rVX.H1'])){
    # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Vaccination effects')
    # text(.5, .5, 'NA', cex = 2)
  }else{
    plot(18:90, colMeans(vaccination.H1), col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'age', main = 'Vaccination effects', ylim = c(.5, 1.5), bty = 'n')
    points(18:90, colMeans(vaccination.H3), col = 'firebrick1', cex = .7)
    abline(h = 1, lty = 2)
  }
  mtext(text = 'C', side = 3, line = 1, at = 10, font = 2)


  ## Imprinting protection
  if(is.na(pars['rPro.H1'])){
    # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Imprinting effects')
    # text(.5, .5, 'NA', cex = 2)
  }else{
    ## Pull out one row for imprinting protection. Focus on Poland, 16.17
    plot(-(2017-18:90), imprinting.H1[58, ], col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'birth year', main = 'Imprinting effects', ylim = c(0,1), xaxt = 'n', bty = 'n')
    axis(side = 1, at = seq(-2017, -1927, by = 10), labels = seq(2017, 1927, by = -10), line = 0)
    points(-(2017-18:90), imprinting.H3[58, ], col = 'firebrick1', cex = .7)
    abline(h = 1, lty = 2)
  }
  mtext(text = 'D', side = 3, line = 1, at = 10, font = 2)


  ## Relative risk
  parnames = c("ages 18-24",  "ages 25-31", "fixed baseline",  "ages 39-45",  "ages 46-52",  "ages 53-59",  "ages 60-66", "ages 67-73",  "ages 74-80",  "ages 81-90", "", "", "", "antivirals", "underlying", "vaccination", "imprinting")
  xvals = c(pars[grep(pattern = "r\\d\\d.\\d\\d", x = names(pars))], NA, NA, NA, pars['rAV'],pars['rDX'], pars['rVX.H1'], pars['rVX.H3'], pars['rPro.H1'], pars['rPro.H3'], NA)
  ## Add CIs
  CIord = sapply(c(grep(pattern = "r\\d\\d.\\d\\d", x = CIplot$parname, value = TRUE), NA, NA, NA,'rAV','rDX','rVX.H1','rVX.H3','rPro.H1', 'rPro.H3', NA), FUN = function(ss) {
    oo = which(CIplot$parname == ss)
    if(length(oo)==0){oo = NA}
    oo
    })
  xlows = CIplot$low[CIord]
  xhis = CIplot$high[CIord]
  yvals = c(1, 2, 4:13, 14, 15, 16-.1, 16+.1, 17-.1, 17+.1, 3)
  par(mar = c(4, 7, 3, 1))
  plot(xvals, yvals, xlim = c(0, 1.6), xaxt = 'n', yaxt = 'n', xlab = 'Relative risk estimate', ylab = '', col = c(rep('black', 14), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), pch = 9, main = 'Parameter estimates')
  segments(x0 = xlows, y0 = yvals, x1 = xhis, col = c(rep('black', 14), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), lwd = 2)
  axis(side = 1, at = seq(0, 1.6, by = .25))
  axis(side = 2, at = 1:17, labels = parnames, las = 2)
  abline(v = 1, lty = 2)
  mtext(text = 'E', side = 3, line = 1, at = -.5, font = 2)


## This old version returns age-specific predictions
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = tested.master/rowSums(tested.master)*age.risk * antivirals.H1 * underlying.H1 * vaccination.H1 * imprinting.H1
  pp.H3 = tested.master/rowSums(tested.master)* age.risk * antivirals.H3 * underlying.H3 * vaccination.H3 * imprinting.H3
  
  ## Normalize,   ## Scale by number of observed cases
  pp.H1 = pp.H1/rowSums(pp.H1)*rowSums(H1.master)
  pp.H3 = pp.H3/rowSums(pp.H3)*rowSums(H3.master)

  
  ##### CONVERT to  BIRTH YEAR
  cur_year = 2000+as.numeric(sub(pattern = '.+(\\d{2}$)', replacement = "\\1", x = rownames(pp.H1)))
  ## Generate a matrix that converts age-specific indices to birth-year sp indices
  bymat = t(sapply(cur_year, function(cc) -(18:90)+cc))
  byind = bymat-1919
  ## Set up empty matrices to store birth year-specific predictions
  H3.converted = H1.converted = H1dat.converted = H3dat.converted = matrix(0, nrow = nrow(pp.H1), ncol = max(byind), dimnames = list(rownames(pp.H1), 1920:1999))
  for(rr in 1:nrow(H1.converted)){
    H1.converted[rr,byind[rr,]] = pp.H1[rr,]
    H3.converted[rr,byind[rr,]] = pp.H3[rr,]
    H1dat.converted[rr, byind[rr,]] = H1.master[rr,]
    H3dat.converted[rr, byind[rr,]] = H3.master[rr,]
  }
  
  # any(rowSums(pp.H1)-rowSums(H1.converted) !=0)
  # any(rowSums(pp.H3)-rowSums(H3.converted)!=0)
  

  return(rbind(colSums(H1.converted), colSums(H3.converted), colSums(H1dat.converted), colSums(H3dat.converted)))
} ##### END FUNCTION
## ------------------



####### GET PREDICTIONS FOR THE FIVE BEST MODELS
par(mfrow = c(3,2))    
pred1 = plotmod(modname = gsub(pattern = 'lk.(\\w+)', '\\1', x = names(del.AIC)[1])); par(mfrow = c(3,2))
pred2 = plotmod(modname = gsub(pattern = 'lk.(\\w+)', '\\1', x = names(del.AIC)[2])); par(mfrow = c(3,2))
pred3 = plotmod(modname = gsub(pattern = 'lk.(\\w+)', '\\1', x = names(del.AIC)[3])); par(mfrow = c(3,2))
pred4 = plotmod(modname = gsub(pattern = 'lk.(\\w+)', '\\1', x = names(del.AIC)[4])); par(mfrow = c(3,2))
pred5 = plotmod(modname = gsub(pattern = 'lk.(\\w+)', '\\1', x = names(del.AIC)[5]))
prednull = plotmod(modname = 'A')






{
  ###########################################################################
  ######################### START INSIGHT PDF ###############################
  ###########################################################################
  ## ------------------------------------------------------------------------
  pdf(outfile2, height = 8)
  #layout(mat = matrix(c(1,1,6,6,6,6, 2,2,6,6,6,6, 3,3,4,4,5,5, 7,7,7,8,8,8), ncol = 6, byrow = T), heights = c(1,1,1,.8))
  layout(mat = matrix(c(1,1,2,2,5,5, 3,3,4,4,5,5, 6,6,6,7,7,7, 6,6,6,8,8,8), ncol = 6, byrow = T), heights = c(1,1,1.4, 1.4))
  par(mgp = c(2,1,0))
  par(mar = c(4,3,3,1))
  #layout.show(6)
  ##########################
  ## 1-3
  ##########################
  npred = plotmod(modname = gsub(pattern = 'lk.(\\w+)', '\\1', x = names(del.AIC)[1])) # Get predicted age distributions 
  
  ##########################
  ## 4. Plot AIC rankings
  ##########################
  par(mar = c(0,5,3,0))
  ## Write a function to plot a filled rectangle of fixed size, given the center coordinate
  rr = function(center, width = 1, height = 1, col.in = 'navy'){
    rect(xleft = center[1]-width/2, ybottom = center[2]-height/2, xright = center[1]+width/2, ytop = center[2]+height/2, col = col.in, border = 'white')
  }
  ## Repeat, only plotting non-degenerate models
  ## Plot AIC results
  mods = fits
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
  for(ii in T.valid){ rr(c(1,ii), col.in = cols1[ii])}
  for(ii in U.valid){ rr(c(2,ii), col.in = cols1[ii])}
  for(ii in V.valid){ rr(c(3,ii), col.in = cols1[ii])}
  for(ii in S.valid){ rr(c(4,ii), col.in = cols1[ii])}
  for(ii in G.valid){ rr(c(5,ii), col.in = cols1[ii])}
  for(ii in N.valid){ rr(c(6,ii), col.in = cols1[ii])}
  for(ii in 1:length(del.AIC)){
    text(7, ii, paste(round(rev(del.AIC)[ii], 2)))
  }
  mtext('F', side = 3, line = 1.5, at = -.25, font = 2)
  
  
  ##########################
  ## 5-6. Plot model fits to observed data
  ##########################
  par(mar = c(3,3,.5,2)+.5)
  cols3 = rev(cols1)
  plot(-1920:-1999, prednull[3,], type = 'l', lwd = 2, xlab = 'birth year', ylab = 'case count', main = 'Model fits: H1N1', xaxt = 'n')
  axis(side = 1, at = seq(-1920, -1999, by = -10), labels = seq(1920, 1999, by = 10), line = 0)
  lines(-1920:-1999, prednull[1,], col = 'pink', lty = 1, lwd = 1.5)
  lines(-1920:-1999, pred4[1,], col = cols3[4], lwd = 1.5)
  lines(-1920:-1999, pred3[1,], col = cols3[3], lwd = 1.5)
  lines(-1920:-1999, pred2[1,], col = cols3[2], lwd = 1.5)
  lines(-1920:-1999, pred1[1,], col = cols3[1], lwd = 1.5)
  legend('topright', legend = c('observed', gsub(pattern = "lk.(\\w+)", replacement = "\\1", names(del.AIC)[1:4]), 'A   '), col = c('black', cols3[1:4], 'pink'), lwd = 1.5, bty = 'n')
  mtext('G', side = 3, line = .5, at = -2000, font = 2)
  
  plot(-1920:-1999, prednull[4,], type = 'l', lwd = 2, xlab = 'birth year', ylab = 'case count', main = 'Model fits: H3N2', xaxt = 'n')
  axis(side = 1, at = seq(-1920, -1999, by = -10), labels = seq(1920, 1999, by = 10), line = 0)
  lines(-1920:-1999, prednull[2,], col = 'pink', lty = 1, lwd = 1)
  lines(-1920:-1999, pred4[2,], col = cols3[4], lwd = 1.5)
  lines(-1920:-1999, pred3[2,], col = cols3[3], lwd = 1.5)
  lines(-1920:-1999, pred2[2,], col = cols3[2], lwd = 1.5)
  lines(-1920:-1999, pred1[2,], col = cols3[1], lwd = 1.5)
  mtext('H', side = 3, line = .5, at = -2000, font = 2)
  dev.off()
  ## ------------------------------------------------------------------------
  ###########################################################################
  ###########################################################################
}












# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### Diagnose INSIGHT fits -- why are different imprinting fits not different?
# plotmod_diag = function(pars, pro.H1, pro.H3){
#   # 1. Assign parameters to be fit
#   rAV = ifelse(is.na(pars['rAV']), 1, pars['rAV']) # Relative risk given antiviral use
#   rDX = ifelse(is.na(pars['rDX']), 1, pars['rDX']) # Relative risk given underlying symptoms
#   rVX.H1 = ifelse(is.na(pars['rVX.H1']), 1, pars['rVX.H1']) # Relative risk given vaccination
#   rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
#   rVX.H3 = ifelse(is.na(pars['rVX.H3']), 1, pars['rVX.H3']) # Relative risk given vaccination
#   rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
#   r18.24 = pars['r18.24'] # Baseline expectation for 18-24 year olds
#   r25.31 = pars['r25.31']
#   b =1 # Fix this as a free paramter. Then estimate all others as relative risk. Most should be lower, bounded at 0.
#   r39.45 = pars['r39.45']
#   r46.52 = pars['r46.52']
#   r53.59 = pars['r53.59']
#   r60.66 = pars['r60.66']
#   r67.73 = pars['r67.73']
#   r74.80 = pars['r74.80']
#   r81.90 = pars['r81.90']
#   
#   ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
#   age.risk = (b*r18.24*a18.24+b*r25.31*a25.31+ b*a32.38+ b*r39.45*a39.45+ b*r46.52*a46.52+ b*r53.59*a53.59+ b*r60.66*a60.66+ b*r67.73*a67.73+ b*r74.80*a74.80+ b*r81.90*a81.90)
#   age.rr = age.risk # save non-normalized version for plot
#   age.risk = age.risk/rowSums(age.risk)
#   
#   antivirals.H1 = (av.master*rAV+(1-av.master))
#   antivirals.H3 = (av.master*rAV+(1-av.master))
#   antivirals = rbind(antivirals.H1, antivirals.H3)
#   
#   underlying.H1 = (dx.master*rDX+(1-dx.master))
#   underlying.H3 = (dx.master*rDX+(1-dx.master))
#   underlying = rbind(underlying.H1, underlying.H3)
#   
#   vaccination.H1 = (vac.master*rVX.H1+(1-vac.master))
#   vaccination.H3 = (vac.master*rVX.H3+(1-vac.master))
#   
#   imprinting.H1 = (pro.H1*rPro.H1+(1-pro.H1))
#   imprinting.H3 = (pro.H3*rPro.H3+(1-pro.H3))
#   
#   
#   ## All rows in age baseline are the same, so plot one arbitrarily
#   plot(18:90, age.rr[1,], main = 'Age effects', ylab = 'relative risk', xlab = 'age', ylim = c(.5, 1.5))
#   abline(h = 1, lty = 2)
#   mtext(text = 'A', side = 3, line = 1, at = 10, font = 2)
#   
#   ## Antivial effects
#   if(is.na(pars['rAV'])){
#     # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Antiviral effects')
#     # text(.5, .5, 'NA', cex = 2)
#   }else{
#     plot(18:90, colMeans(antivirals), col = 'black', cex = .7, ylab = 'relative risk', xlab = 'age', main = 'Antiviral effects', ylim = c(.5, 1.5))
#     abline(h = 1, lty = 2)
#   }
#   mtext(text = 'B', side = 3, line = 1, at = 10, font = 2)
#   
#   ## Underlying symptoms effects
#   if(is.na(pars['rDX'])){
#     # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Underlying cond. effects')
#     # text(.5, .5, 'NA', cex = 2)
#   }else{
#     plot(18:90, colMeans(underlying), col = 'black', cex = .7, ylab = 'relative risk', xlab = 'age', main = 'Underlying cond. effects', ylim = c(.5, 1.5))
#     abline(h = 1, lty = 2)
#     mtext(text = 'C', side = 3, line = 1, at = 10, font = 2)
#   }
#   
#   
#   ## Vaccination
#   if(is.na(pars['rVX.H1'])){
#     # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Vaccination effects')
#     # text(.5, .5, 'NA', cex = 2)
#   }else{
#     plot(18:90, colMeans(vaccination.H1), col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'age', main = 'Vaccination effects', ylim = c(.5, 1.5))
#     points(18:90, colMeans(vaccination.H3), col = 'firebrick1', cex = .7)
#     abline(h = 1, lty = 2)
#   }
#   mtext(text = 'C', side = 3, line = 1, at = 10, font = 2)
#   
#   
#   ## Imprinting protection
#   if(is.na(pars['rPro.H1'])){
#     # plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = 'Imprinting effects')
#     # text(.5, .5, 'NA', cex = 2)
#   }else{
#     ## Pull out one row for imprinting protection. Focus on Poland, 16.17
#     plot(2017-18:90, imprinting.H1[58, ], col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'birth year', main = 'Imprinting effects', ylim = c(0,1))
#     points(2017-18:90, imprinting.H3[58, ], col = 'firebrick1', cex = .7)
#     abline(h = 1, lty = 2)
#   }
#   mtext(text = 'D', side = 3, line = 1, at = 10, font = 2)
#   
#   
#   ## Relative risk
#   parnames = c("ages 18-24",  "ages 25-31", "fixed baseline",  "ages 39-45",  "ages 46-52",  "ages 53-59",  "ages 60-66", "ages 67-73",  "ages 74-80",  "ages 81-90", "", "", "", "antivirals", "underlying", "vaccination", "imprinting")
#   xvals = c(pars[grep(pattern = "r\\d\\d.\\d\\d", x = names(pars))], NA, NA, NA, pars['rAV'],pars['rDX'], pars['rVX.H1'], pars['rVX.H3'], pars['rPro.H1'], pars['rPro.H3'], NA)
#   ## Add CIs
#   # CIord = sapply(c(grep(pattern = "r\\d\\d.\\d\\d", x = CIplot$parname, value = TRUE), NA, NA, NA,'rAV','rDX','rVX.H1','rVX.H3','rPro.H1', 'rPro.H3', NA), FUN = function(ss) {
#   #   oo = which(CIplot$parname == ss)
#   #   if(length(oo)==0){oo = NA}
#   #   oo
#   # })
#   # xlows = CIplot$low[CIord]
#   # xhis = CIplot$high[CIord]
#   yvals = c(1, 2, 4:13, 14, 15, 16-.1, 16+.1, 17-.1, 17+.1, 3)
#   par(mar = c(4, 7, 3, 1))
#   plot(xvals, yvals, xlim = c(0, 1.4), xaxt = 'n', yaxt = 'n', xlab = 'Relative risk estimate', ylab = '', col = c(rep('black', 14), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), pch = 9, main = 'Parameter estimates')
#  # segments(x0 = xlows, y0 = yvals, x1 = xhis, col = c(rep('black', 14), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), lwd = 2)
#   axis(side = 1, at = seq(0, 1.25, by = .25))
#   axis(side = 2, at = 1:17, labels = parnames, las = 2)
#   abline(v = 1, lty = 2)
#   mtext(text = 'E', side = 3, line = 1, at = -.5, font = 2)
#   
#   
#   ## This old version returns age-specific predictions
#   # 2. calculate predicted distribution, pp, as a function of the parameters:
#   # This step gives the model prediction
#   pp.H1 = tested.master/rowSums(tested.master)*age.risk * antivirals.H1 * underlying.H1 * vaccination.H1 * imprinting.H1
#   pp.H3 = tested.master/rowSums(tested.master)* age.risk * antivirals.H3 * underlying.H3 * vaccination.H3 * imprinting.H3
#   
#   ## Normalize,   ## Scale by number of observed cases
#   pp.H1 = pp.H1/rowSums(pp.H1)*rowSums(H1.master)
#   pp.H3 = pp.H3/rowSums(pp.H3)*rowSums(H3.master)
#   
#   
#   ##### CONVERT to  BIRTH YEAR
#   cur_year = 2000+as.numeric(sub(pattern = '.+(\\d{2}$)', replacement = "\\1", x = rownames(pp.H1)))
#   ## Generate a matrix that converts age-specific indices to birth-year sp indices
#   bymat = t(sapply(cur_year, function(cc) -(18:90)+cc))
#   byind = bymat-1919
#   ## Set up empty matrices to store birth year-specific predictions
#   H3.converted = H1.converted = H1dat.converted = H3dat.converted = matrix(0, nrow = nrow(pp.H1), ncol = max(byind), dimnames = list(rownames(pp.H1), 1920:1999))
#   for(rr in 1:nrow(H1.converted)){
#     H1.converted[rr,byind[rr,]] = pp.H1[rr,]
#     H3.converted[rr,byind[rr,]] = pp.H3[rr,]
#     H1dat.converted[rr, byind[rr,]] = H1.master[rr,]
#     H3dat.converted[rr, byind[rr,]] = H3.master[rr,]
#   }
#   
#   # any(rowSums(pp.H1)-rowSums(H1.converted) !=0)
#   # any(rowSums(pp.H3)-rowSums(H3.converted)!=0)
#   
#   
#   return(rbind(colSums(H1.converted), colSums(H3.converted), colSums(H1dat.converted), colSums(H3dat.converted)))
# } ##### END FUNCTION
# ## ------------------
# 
# 
# best = plotmod_diag(pars = fits$lk.ATVN$par, pro.H1 = proN1.master, pro.H3 = proN2.master)
# ## Change the protection input, but not the par values and see what happend
# best_S = plotmod_diag(pars = fits$lk.ATVN$par, pro.H1 = proH1.master, pro.H3 = proH3.master)
# best_G = plotmod_diag(pars = fits$lk.ATVN$par, pro.H1 = prog1.master, pro.H3 = prog2.master)
# best_A = plotmod_diag(pars = fits$lk.ATVN$par, pro.H1 = prog1.master*0, pro.H3 = prog2.master*0)
# par(mfrow = c(2,1))
# cols = viridis_pal(alpha = .9)(4)
# xx = barplot(best[3,], col = 'gray80', border = NA, space = 0, main = 'H1N1', ylim = c(0, 70))
# lines(xx, best[1,], col = cols[1])
# lines(xx, best_S[1,], col = cols[2])
# lines(xx, best_G[1,], col = cols[3])
# lines(xx, best_A[1,], col = cols[4])
# 
# xx = barplot(best[4,], col = 'gray80', border = NA, space = 0, main = 'H3N2', ylim = c(0, 70))
# lines(xx, best[2,], col = cols[1])
# lines(xx, best_S[2,], col = cols[2])
# lines(xx, best_G[2,], col = cols[3])
# lines(xx, best_A[2,], col = cols[4])
# legend('topleft',c('best', 'S', 'G', 'no imprinting'), col = cols, lty = 1, bty = 'n')
# 
# 
# 
# 
# b2nd = plotmod_diag(pars = fits$lk.ATVS$par, pro.H1 = proH1.master, pro.H3 = proH3.master)
# AA = plotmod_diag(pars = fits$lk.ATV$par, pro.H1 = 1, pro.H3 = 1)
# par(mfrow = c(1,1))
# cols = viridis_pal(alpha = .9)(4)
# xx = barplot(best[3,], col = 'gray80', border = NA, space = 0, main = 'H1N1', ylim = c(0, 70))
# lines(xx, best[1,], col = cols[1])
# lines(xx, b2nd[1,], col = cols[2])
# lines(xx, AA[1,], col = cols[3])
# 
# xx = barplot(best[4,], col = 'gray80', border = NA, space = 0, main = 'H3N2', ylim = c(0, 70))
# lines(xx, best[2,], col = cols[1])
# lines(xx, b2nd[2,], col = cols[2])
# lines(xx, AA[2,], col = cols[3])
# legend('topleft',c('best', '2nd best'), col = cols[1:2], lty = 1, bty = 'n')
# 
