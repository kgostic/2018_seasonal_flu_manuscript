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
setwd('~/Dropbox/R/2018_seasonal_flu/2017_AZ/')
load('processed-data/AZ_model_fits.RData')
load('processed-data/AZ_CIs.RData')
source('00-Inputs_multinomial.R') 
setwd('../')# Switch back into the home directory


## OUTPUTS
plot3 = 'figures/AZ_H1N1_fit.pdf'
plot4 = 'figures/AZ_H3N2_fit.pdf'
plot5 = 'figures/AZ_age_fit.pdf'
plot6 = 'figures/AZ_imp_fit.pdf'




## Set up color palette
cols = rev(viridis_pal(alpha = .7, begin = 1, end = .1)(5))
show_col(cols[1])
cols[2] = 'purple'
cols[3] = 'dodgerblue'
cols[4] = 'limegreen'
cols[5] = 'goldenrod'

cols = c(NA, 'goldenrod', 'limegreen', 'dodgerblue', 'purple')

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
  plot(0:97, age.baseline.rr[13,], main = 'Age effects', ylab = 'relative risk', xlab = 'age', bty = 'n', ylim = c(0,1))
  mtext('A', side = 3, line = 1.5, at = -5, font = 2)
  #abline(h = 1, lty = 2)
  
  ## Predicted effects from imprinting protection
  ## Plot the colmeans (mean predicted protection across birth years)
  ## This is necessary because imprinting reconstructions differ slightly from year to year, as children (with recent birth years) get older.
  if(is.na(pars['rPro.H1'])){ # If no imprinting protection par, this factor was not relevant to the model fit, and we plot an empty windwo
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste('Imprinting protection\n', i.type, sep = ''), bty = 'n')
    text(.5, .5, 'NA', cex = 2)
  }else{ # Else, if imprinting protection was included in the model, plot mean birth year-specific relative risk, after adjusting for imprinting protection.
    ## Plot just one row here, because predicted protection changes with birth year
    plot(-2015:-1918, imprinting.H1[13,], col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'birth year', main = 'Imprinting effects', ylim = c(0,1), bty = 'n', xaxt = 'n')
    axis(side = 1, at = seq(-2015, -1918, by = 10), labels = seq(2015, 1918, by = -10), line = 0)
    points(-2015:-1918, imprinting.H3[13,], col = 'firebrick1', cex = .7)
    }
  ## H1N1 risk in blue, H3N2 risk in red
  mtext('B', side = 3, line = 1.5, at = -5, font = 2)
  #abline(h = 1, lty = 2)
  
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
  return(list(age = age.baseline.rr, iH1 = imprinting.H1, iH3 = imprinting.H3,
              fits = rbind(rev(colSums(pp.H1/rowSums(pp.H1)*rowSums(H1.master))), rev(colSums(pp.H3/rowSums(pp.H3)*rowSums(H3.master))))))
}
# Get prediction for model AS, but don't save to pdf
AS = plotmod1(lk.AS$par, pro.H1 = proH1.master, pro.H3 = proH3.master, i.type = 'HA Sub', CIs = AS.CIs)
spred = AS$fits# Get predicted age distributions using hte best fit subtype-specific imprinting model.
AG = plotmod1(lk.AG$par, pro.H1 = prog1.master, pro.H3 = prog2.master, i.type = "HA grp", CIs = AG.CIs)
gpred = AG$fits
AA = plotmod1(lk.A$par, pro.H1 = 1, pro.H3 = 1, i.type = NULL, CIs = A.CIs)
apred = AA$fits
AN = plotmod1(lk.AN$par, pro.H1 = proN1.master, pro.H3 = proN2.master, i.type = 'NA Sub', CIs = AN.CIs) 
npred = AN$fits





######### Plot AZ age fits
pdf(file = plot5, width = 3, height = 1.5)
par(mar = c(3,3,1,1)+.5, mgp = c(2,1,0))
plot(0:97, AN$age[13,], main = '', ylab = 'relative risk', xlab = 'age', bty = 'n', ylim = c(0,1))
#abline(h = 1, lty = 2)
dev.off()


######### Plot AZ imprinting fits
pdf(file = plot6, width = 3, height = 1.5)
par(mar = c(3,3,1,1)+.5, mgp = c(2,1,0))
plot(-2015:-1918, AN$iH1[13,], col = 'dodgerblue', cex = .7, ylab = 'relative risk', xlab = 'birth year', main = '', ylim = c(0,1), bty = 'n', xaxt = 'n')
axis(side = 1, at = seq(-2015, -1918, by = 10), labels = seq(2015, 1918, by = -10), line = 0)
points(-2015:-1918, AN$iH3[13,], col = 'firebrick1', cex = .7)
dev.off()


######### Plot AZ H1N1 fits
dal = round(del.AIC, 2); names(dal) = gsub(pattern = 'lk.(\\w+)', replacement = '\\1', x = names(del.AIC))
pdf(file = plot3, width = 4, height = 3.5)
par(mar = c(3,3,2,2)+.5, mgp = c(2,1,0))
xx = barplot(colSums(H1.master), col = 'gray', border = 'gray', ylim = c(0, 115), xlab = 'birth year', ylab = 'cases')
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA, line = 0)
lines(xx, rev(apred[1,]), col = cols[2], lwd = 1.7)
lines(xx, rev(gpred[1,]), col = cols[3], lwd = 1.7)
lines(xx, rev(spred[1,]), col = cols[4], lwd = 1.7)
lines(xx, rev(npred[1,]), col = cols[5], lwd = 1.7)
legend('topright', legend = c('observed', expression(paste('AN fit, ', Delta, 'AIC=', 0.00)), 
                              expression(paste('AS fit, ', Delta, 'AIC=', 23.42)),
                              expression(paste('AG fit, ', Delta, 'AIC=', 245.18)),
                              expression(paste('A   fit, ', Delta, 'AIC=', 380.47))), pch = c(15, NA, NA, NA, NA), col = c('gray', cols[5:2]), lwd = c(NA,2,2,1,1), bty = 'n')
dev.off()


######### Plot AZ H3N2 fits
pdf(file = plot4, width = 4, height = 3.5)
par(mar = c(3,3,2,2)+.5, mgp = c(2,1,0))
xx = barplot(colSums(H3.master), col = 'gray', border = 'gray', xlab = 'birth year', ylab = 'cases')
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA, line = 0)
lines(xx, rev(apred[2,]), col = cols[2], lwd = 1.7)
lines(xx, rev(gpred[2,]), col = cols[3], lwd = 1.7)
lines(xx, rev(spred[2,]), col = cols[4], lwd = 1.7)
lines(xx, rev(npred[2,]), col = cols[5], lwd = 1.7)
legend('topright', legend = c('observed', expression(paste('AN fit, ', Delta, 'AIC=', 0.00)), 
                              expression(paste('AS fit, ', Delta, 'AIC=', 23.42)),
                              expression(paste('AG fit, ', Delta, 'AIC=', 245.18)),
                              expression(paste('A   fit, ', Delta, 'AIC=', 380.47))), pch = c(15, NA, NA, NA, NA), col = c('gray', cols[5:2]), lwd = c(NA,2,2,1,1), bty = 'n')
dev.off()


