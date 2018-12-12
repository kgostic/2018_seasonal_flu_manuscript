## Import data and run multinomila model comparison
## Fit to AZ data, from all seasons
## Fit a single age curve to all data
## Test three impriting hypothese: HA group level, HA subtype level, NA subtype level
## PLOT RESULTS
## SAVE MODEL FITS 

## Clear memory
rm(list = ls())
#setwd('2017_AZ/')


## OUTPUTS
modelfits = 'processed-data/AZ_model_fits.RData'
outfile1 = '../figures/AZ_predictions.pdf' # Plot of model predictions vs. observed data
outfile2 = '../figures/AZ_AIC.pdf' # Plot of AIC scores and factors included
outfile3 = '../figures/AZ_NA_model_results.pdf' # Plot of AIC scores and factors included
outfile4 = '../figures/AZ_HAsub_model_results.pdf' # Plot of AIC scores and factors included

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





## Define a function analagous to the likelihood above, which outputs model predictions
#### INPUTS
##     pars - named vector of paramter values (should be best estimates from fits above)
##     pro.H1 - probabilities of H1N1 protection
##     pro.H3 - probabilities of H3N2 protection
##     i.type - character, type of imprinting protection (can be 'HA subtype', 'HA group', or 'NA subtype'). Used in plot titles.
plotmod1 = function(pars, pro.H1 = 1, pro.H3 = 1, i.type = NULL){
  
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
  age.baseline = age.baseline/rowSums(age.baseline)
  imprinting.H1 = (pro.H1*rPro.H1+(1-pro.H1))
  imprinting.H3 = (pro.H3*rPro.H3+(1-pro.H3))
  
  
  ## Plot predictions vs. observed data
  par(mfrow = c(2,2))
  
  ## Age-specific risk prediction
  ## Plot row 15, which represents the last observed season, and aligns 0-year-olds with the first column (2015 birth year)
  ## No need to take the mean, because age-specific predictions are identical across countries and years
  plot(0:97, age.baseline[15,], main = 'Age effects', ylab = 'proportion of cases', xlab = 'case age')
  
  
  ## Predicted effects from imprinting protection
  ## Plot the colmeans (mean predicted protection across birth years)
  ## This is necessary because imprinting reconstructions differ slightly from year to year, as children (with recent birth years) get older.
  if(is.na(pars['rPro.H1'])){ # If no imprinting protection par, this factor was not relevant to the model fit, and we plot an empty windwo
    plot(1, 1, xlim = c(0,1), ylim = c(0,1), col = 'white', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '', main = paste('Imprinting protection\n', i.type, sep = ''))
    text(.5, .5, 'NA', cex = 2)
  }else{ # Else, if imprinting protection was included in the model, plot mean birth year-specific relative risk, after adjusting for imprinting protection.
    plot(colMeans(imprinting.H1), col = 'dodgerblue', cex = .7, ylab = 'est. relative risk impact', xlab = 'case birth year', main = 'Imprinting', ylim = c(0,1), xaxt = 'n')
    points(colMeans(imprinting.H3), col = 'firebrick1', cex = .7)
    axis(1, at = seq(3, 98, by = 10), labels = seq(1920, 2015, by = 10), las = 2) }
  ## H1N1 risk in blue, H3N2 risk in red
  
  
  ## Plot relative risk point estimates for each free parameter.
  parnames = c("ages 0-4", "ages5-10", "ages 11-17", "ages 18-24",  "ages 25-31", "ages 32-38",  "ages 39-45",  "ages 46-52",  "ages 53-59",  "ages 60-66", "ages 67-73",  "ages 74-80",  "ages 81+", "", "", "", paste('imprinting', i.type, sep = ', ')) 
  xvals = c(pars[grep(pattern = "r\\d\\d?.\\d\\d", x = names(pars))], NA, NA, NA, pars['rPro.H1'], pars['rPro.H3'], 1)
  yvals = c(1:15, 16, 16, 0)+1
  par(mar = c(4, 7, 2, 1))
  plot(xvals, yvals, xlim = c(0, 1.3), xaxt = 'n', yaxt = 'n', xlab = 'Relative risk estimate', ylab = '', col = c(rep('black', 12), 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'dodgerblue', 'firebrick1', 'black'), pch = 9)
  axis(side = 1, at = seq(0, 1.25, by = .25))
  axis(side = 2, at = 1:17, labels = parnames, las = 2)
  abline(v = 1, lty = 2)
  
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction
  pp.H1 = age.baseline * imprinting.H1
  pp.H3 = age.baseline * imprinting.H3
  ## Return the predicted age distributions of infection. Plot these below against the observed data
  return(rbind(colSums(pp.H1/rowSums(pp.H1)*rowSums(H1.master)), colSums(pp.H3/rowSums(pp.H3)*rowSums(H3.master))))
}


## Plot best model
par(mfrow = c(4,4,2,3))
pdf(outfile3)
npred = plotmod1(lk.AN$par, pro.H1 = proN1.master, pro.H3 = proN2.master, i.type = 'NA Subtype') # Get predicted age distributions using the best fit NA imprinting model
dev.off()

pdf(outfile4)
spred = plotmod1(lk.AS$par, pro.H1 = proH1.master, pro.H3 = proH3.master, i.type = 'HA Subtype') # Get predicted age distributions using hte best fit subtype-specific imprinting model.
dev.off()

## Plot model fits to data
pdf(outfile1, width = 7, height = 4.5)
par(mfrow = c(1,2))
plot(2015:1918, colSums(H1.master), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H1N1')
lines(2015:1918, npred[1,], col = 'magenta', lwd = 2)
lines(2015:1918, spred[1,], col = 'orange', lwd = 2)
legend('topleft', legend = c('N', 'S'), col = c('magenta', 'orange'), lty = 1)

plot(2015:1918, colSums(H3.master), type = 'l', lwd = 2, xlab = 'birth year', ylab = 'predicted total case count', main = 'H3N2')
lines(2015:1918, npred[2,], col = 'magenta', lwd = 2)
lines(2015:1918, spred[2,], col = 'orange', lwd = 2)
dev.off()




## Plot ranked model fits, in terms of AIC value
pdf(outfile2, width = 3.25, height = 2)
par(mar = c(2,3,2,2))
## Write a function to plot a filled rectangle of fixed size, given the center coordinate
rr = function(center, width = 1, height = 1, col.in = 'navy'){
  rect(xleft = center[1]-width/2, ybottom = center[2]-height/2, xright = center[1]+width/2, ytop = center[2]+height/2, col = col.in, border = 'black')
}
## Repeat, only plotting non-degenerate models
## Plot AIC results
mods2 = mods[rev(names(del.AIC))]
factors = c('A', 'S', 'G', 'N')
S.valid = grep('S', names(mods2))
G.valid = grep('G', names(mods2))
N.valid = grep('N', names(mods2))
A.valid = 1:length(mods2)
plot.new()
plot.window(xlim = c(0.5, 5.5), ylim = c(0.5, length(del.AIC)+.5))
axis(2, at = 1:length(mods2), labels = gsub(pattern = 'lk.(\\w+)',replacement = "\\1", names(mods2)), las = 2)
axis(3, at = 1:5, labels = c(factors, expression(paste(Delta, 'AIC', sep = ''))))
for(ii in A.valid){ rr(c(1,ii))}
for(ii in S.valid){ rr(c(2,ii))}
for(ii in G.valid){ rr(c(3,ii))}
for(ii in N.valid){ rr(c(4,ii))}
for(ii in 1:length(del.AIC)){
text(5, ii, paste(round(rev(del.AIC)[ii], 2)))
}
dev.off()



