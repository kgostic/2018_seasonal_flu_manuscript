## plot pandemic fits

##########################################
## Start with results from AZ data
##########################################
rm(list = ls())

outfile1 = 'figures/AZ_pandemic_fits.pdf'
outfile2 = 'figures/AZ_pandemic_resid.pdf'

## Temporarily switch to the 2017_AZ subdirectory and load AZ model fits, CIs and model inputs
## We will load the INISHGT inputs below, when we start that plot.
setwd('2017_AZ/')
load('processed-data/AZ_model_fits.RData')
load('processed-data/AZ_pandemic_fits.RData')
load('processed-data/AZ_CIs.RData')
load('processed-data/AZ_profiles.RData')
source('00-Inputs_multinomial.R') 
setwd('../')# Switch back into the home directory

pdf(outfile1, width = 3)
par(mfrow = c(3,1), mgp = c(2,1,0), mar = c(4,3,3,1))
## Plot the estimated relative risk given imprinting, as fitted to seasonal data
plot(AZ_seasonal_ests, 1:3, col = c('limegreen', 'purple', 'dodgerblue'), pch = 16, xlim = c(0, 1), ylim = c(.5, 3.5), xaxt = 'n', yaxt = 'n', ylab = '', bty = 'n', xlab = 'Relative risk given imprinting')
axis(1, at = seq(0, 1, by = .2))
axis(2, at = 1:3, labels = c('AG', 'AN', 'AS'))     
points(AZ_firstwave_ests, 1:3, col = c('limegreen', 'purple', 'dodgerblue'), pch = 4)
points(AZ_secondwave_ests, 1:3, col = c('limegreen', 'purple', 'dodgerblue'), pch = 5)

## Add seasonal CIs
segments(x0 = AG.CIs[1,1], y0 = 1, x1 = AG.CIs[2,1], col = 'limegreen')
segments(x0 = AN.CIs[1,1], y0 = 2, x1 = AN.CIs[2,1], col = 'purple')
segments(x0 = AS.CIs[1,1], y0 = 3, x1 = AS.CIs[2,1], col = 'dodgerblue')
legend('topright', c('fitted to seasonal data', 'first pandemic wave', 'second pandemic wave'),  pch = c(16, 4, 5), bty = 'n')
mtext(text = 'A', side = 3, line = 0, at =0, font = 2)


## Define a function to predict the age distribution of cases, given model paramters:
AZ_pdm_prediction = function(pars, fitted.age.pars,wPro.H1){
  dat.H1 = H1.master_2009
  a0.4 = a0.4_2009
  a5.10 = a5.10_2009
  a11.17 = a11.17_2009
  a18.24 = a18.24_2009
  a25.31 = a25.31_2009
  a32.38 = a32.38_2009
  a39.45 = a39.45_2009
  a46.52 = a46.52_2009
  a53.59 = a53.59_2009
  a60.66 = a60.66_2009
  a67.73 = a67.73_2009
  a74.80 = a74.80_2009
  a81.90plus = a81.90plus_2009
  # 1. Assign parameters to be fit (all age paramters, and those named in pars)
  rPro.H1 = (pars['rPro.H1'])# Relative risk given imprinting protection
  b = 1 # Fix relative risk in the baseline group (Ages 0-4) at value 1. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r5.10 = fitted.age.pars['r5.10'] # Relative risk for 5 to 10 year olds (free paramter to estiamte)
  r11.17 = fitted.age.pars['r11.17'] # Relative risk for 11-17 year olds
  r18.24 = fitted.age.pars['r18.24'] # etc.
  r25.31 = fitted.age.pars['r25.31']
  r32.38 = fitted.age.pars['r32.38']
  r39.45 = fitted.age.pars['r39.45']
  r46.52 = fitted.age.pars['r46.52'] 
  r53.59 = fitted.age.pars['r53.59']
  r60.66 = fitted.age.pars['r60.66'] 
  r67.73 = fitted.age.pars['r67.73'] 
  r74.80 = fitted.age.pars['r74.80'] 
  r81.90p = fitted.age.pars['r81.90p'] 
  
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.baseline = b*(a0.4 +
                      r5.10*a5.10 +
                      r11.17*a11.17+
                      r18.24*a18.24+
                      r25.31*a25.31+
                      r32.38*a32.38+ 
                      r39.45*a39.45+ 
                      r46.52*a46.52+ 
                      r53.59*a53.59+ 
                      r60.66*a60.66+ 
                      r67.73*a67.73+ 
                      r74.80*a74.80+ 
                      r81.90p*a81.90plus)
  if(is.null(dim(age.baseline))){ # If only one row
    age.baseline = age.baseline/sum(age.baseline)
  }else{
    age.baseline = age.baseline/rowSums(age.baseline)
  }
  
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction for H1N1 cases
 ppH1 = age.baseline * (wPro.H1*rPro.H1+(1-wPro.H1))
 ppH1 = (ppH1/rowSums(ppH1)*rowSums(H1.master_2009))
 ppH1 = ppH1[,as.character(2009:1918)]
 return(ppH1)
}



## Plot the fits to data
xx = barplot(H1.master_2009[1,as.character(2009:1918)], xlab = 'birth year', ylab = 'cases', col = 'gray', border = 'gray')
lines(xx, AZ_pdm_prediction(pars = firstwave_AS$par, fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009)[1,], col = 'dodgerblue')
lines(xx, AZ_pdm_prediction(pars = firstwave_AN$par, fitted.age.pars = lk.AN$par, wPro.H1 = proN1.master_2009)[1,], col = 'purple')
#lines(2009:1918, AZ_pdm_prediction(pars = lk.AS$par['rPro.H1'], fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009)[1,], col = 'brown1')
lines(xx, AZ_pdm_prediction(pars = firstwave_AG$par, fitted.age.pars = lk.AG$par, wPro.H1 = prog1.master_2009)[1,], col = 'limegreen')
#lines(xx, AZ_pdm_prediction(pars = firstwave_AG$par*0+1, fitted.age.pars = lk.A$par, wPro.H1 = prog1.master_2009)[1,], col = 'yellow')
legend('topright', c('data', expression(paste('AG fit, ', Delta, 'AIC=0.00')),
                     expression(paste('AS fit, ', Delta, 'AIC=20.51')),
                     expression(paste('AN fit, ', Delta, 'AIC=41.87'))),
       col = c('gray', 'limegreen','dodgerblue', 'purple'), lwd = c(NA,1,1,1), bty = 'n', pch = c(15, NA, NA, NA))
mtext(text = 'B', side = 3, line = 0, at = xx[1]-3, font = 2)



xx = barplot(H1.master_2009[2,as.character(2009:1918)], xlab = 'birth year', ylab = 'cases', col = 'gray', border = 'gray')
lines(xx, AZ_pdm_prediction(pars = secondwave_AS$par, fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009)[2,], col = 'dodgerblue')
lines(xx, AZ_pdm_prediction(pars = secondwave_AN$par, fitted.age.pars = lk.AN$par, wPro.H1 = proN1.master_2009)[2,], col = 'purple')
#lines(2009:1918, AZ_pdm_prediction(pars = lk.AS$par['rPro.H1'], fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009)[1,], col = 'brown1')
lines(xx, AZ_pdm_prediction(pars = secondwave_AG$par, fitted.age.pars = lk.AG$par, wPro.H1 = prog1.master_2009)[2,], col = 'limegreen')
#lines(2009:1918, AZ_pdm_prediction(pars = lk.AG$par['rPro.H1'], fitted.age.pars = lk.AS$par, wPro.H1 = proH1.master_2009)[1,], col = 'blue')
legend('topright', c('data', expression(paste('AN fit, ', Delta, 'AIC=0.00')),
                     expression(paste('AS fit, ', Delta, 'AIC=17.69')),
                     expression(paste('AG fit, ', Delta, 'AIC=271.68'))),
                     col = c('gray', 'purple','dodgerblue', 'limegreen'), lwd = c(NA,1,1,1), bty = 'n', pch = c(15, NA, NA, NA))
mtext(text = 'C', side = 3, line = 0, at = xx[1]-3, font = 2)
dev.off()



### Do plot for EEID
{
pdf('EEID_2019.pdf')
par(mfrow = c(3,1), mar = c(3,3,1,1))
xx= barplot(rbind(proH3.master_2009[1,as.character(2009:1918)], 
              proH1.master_2009[1,as.character(2009:1918)],
              prog1.master_2009[1,as.character(2009:1918)]-proH1.master_2009[1,as.character(2009:1918)]),
              col = (c('firebrick1','dodgerblue','lightblue2')), border = NA, space = 0)
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA)
xx= barplot(H1.master_2009[1,as.character(2009:1918)], xlab = 'birth year', ylab = 'cases', col = 'dodgerblue', border = 'dodgerblue')
#lines(xx, AZ_pdm_prediction(pars = firstwave_AS$par*0+1, fitted.age.pars = lk.A$par, wPro.H1 = proH1.master_2009*0+1)[1,], col = 'black')
lines(smooth.spline(xx, colSums(H1.master[,as.character(2009:1918)])/sum(H1.master[,as.character(2009:1918)])*sum(H1.master_2009[1,as.character(2009:1918)])), col = 'black')
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA)

xx= barplot(H1.master_2009[2,as.character(2009:1918)], xlab = 'birth year', ylab = 'cases', col = 'dodgerblue', border = 'dodgerblue')
#lines(xx, AZ_pdm_prediction(pars = firstwave_AS$par*0+1, fitted.age.pars = lk.A$par, wPro.H1 = proH1.master_2009*0+1)[1,], col = 'black')
lines(smooth.spline(xx, colSums(H1.master[,as.character(2009:1918)])/sum(H1.master[,as.character(2009:1918)])*sum(H1.master_2009[2,as.character(2009:1918)])), col = 'black')
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA)
dev.off()
}


## Write the del.AIC values to tables
write.csv(x = round(del.AIC,2), file = 'AZ_seasonal_AICs.csv')
write.csv(x = round(pdm.del.AIC,2), file = 'AZ_pandemic_AICs.csv')



pdf(outfile2, width = 5, height = 5)
par(mfrow = c(2,1), mar = c(3,3,2,2), mgp = c(2,1,0))
AA = AZ_pdm_prediction(pars = firstwave_AG$par*0+1, fitted.age.pars = lk.A$par, wPro.H1 = prog1.master_2009)
xx = barplot(H1.master_2009[1,as.character(2009:1918)]-AA[1,], xlab = 'birth year', ylab = 'excess cases', col = 'gray', border = 'black', ylim = c(-90, 130), main = '')
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA)
pdmxx = xx[which((2009:1918)%in%c(1977, 1968, 1957))]
abline(v = pdmxx, lty = 2)
text(x = pdmxx-2, 100, labels = c('1977', '1968', '1957'), srt = 90)

xx = barplot(H1.master_2009[2,as.character(2009:1918)]-AA[2,], xlab = 'birth year', ylab = 'excess cases', col = 'gray', border = 'black', ylim = c(-90, 130))
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA)
abline(v = pdmxx, lty = 2)
text(x = pdmxx-2, 100, labels = c('1977', '1968', '1957'), srt = 90)
dev.off()

AA = AZ_pdm_prediction(pars = firstwave_AG$par, fitted.age.pars = lk.AG$par, wPro.H1 = prog1.master_2009)[1,]
xx = barplot(H1.master_2009[1,as.character(2009:1918)]-AA, xlab = 'birth year', ylab = 'excess cases', col = 'gray', border = 'black', ylim = c(-90, 130), main = '')
axis(side = 1, at = xx[seq(1, length(xx), by = 10)], labels = NA)
pdmxx = xx[which((2009:1918)%in%c(1977, 1968, 1957))]
abline(v = pdmxx, lty = 2)
text(x = pdmxx-2, 100, labels = c('1977', '1968', '1957'), srt = 90)







##########################################
## Now, do the INSIGHT results
##########################################
rm(list = ls())
## Temporarily switch to the 2017_AZ subdirectory and load AZ model fits, CIs and model inputs
## We will load the INISHGT inputs below, when we start that plot.
setwd('2017_INSIGHT/')
load('processed-data/INSIGHT_fitted_models.RData')
load('processed-data/INSIGHT_fitted_pandemic_models.RData')
load('processed-data/INSIGHT_pandemic_and_seasonal_imp_risk_ests.RData')
CIs = read.csv('processed-data/INSIGHT_par_ests_CIs.csv')
source('00-Import_FLU002_-for-multinomial.R') 
setwd('../')# Switch back into the home directory
par(mfrow = c(1,1))


## extract model names
par(mar = c(3,5,4,1))
modnames = gsub(pattern = 'lk.(\\w+).rPro.H1', replacement = '\\1', x = names(seasonal_ests))
## Re-order by imprinting breadth
ord = c(grep(modnames, pattern = 'G'), grep(modnames, pattern = 'S'), grep(modnames, pattern = 'N'))
modnames = modnames[ord]
seasonal_ests = seasonal_ests[ord]
pandemic_ests = pandemic_ests[ord]
## Plot seasonal estimates
plot(seasonal_ests, 1:length(seasonal_ests), xlim = c(.5, 1), xlab = 'relative risk given imprinting protection', yaxt = 'n', bty = 'n', ylab = '', pch = 16, cex = .7)
axis(side = 2, at = 1:length(seasonal_ests), labels = modnames, las = 2)
## Re-plot the estimates from models containns HA subtype imprinting
sr = grep(modnames, pattern = 'S')
points(seasonal_ests[sr], sr, col = 'orange', pch = 16)
points(pandemic_ests[sr], sr, col = 'orange', pch = 4)
## Repeat from models containns NA subtype imprinting
nr = grep(modnames, pattern = 'N')
points(seasonal_ests[nr], nr, col = 'magenta', pch = 16)
points(pandemic_ests[nr], nr, col = 'magenta', pch = 4)
## Repeat from models containns HA group  imprinting
gr = grep(modnames, pattern = 'G')
points(seasonal_ests[gr], gr, col = 'dodgerblue', pch = 16)
points(pandemic_ests[gr], gr, col = 'dodgerblue', pch = 4)
### Add seasonal CIs
cols = rep(c('dodgerblue', 'orange', 'magenta'), each = 8)
for(ii in 1:length(modnames)){
  CI = as.numeric(CIs[CIs$modname == modnames[ii] & CIs$variable == 'rPro.H1', c('low', 'high')])
  segments(x0 = CI[1], y0 = ii, x1 = CI[2], col = cols[ii])
}
legend(.4, 26.5, c('fitted to seasonal data', 'fitted to pandemic data'),  pch = c(16, 4), bty = 'n', xpd = NA, xjust = 0)



## Write the del.AIC values to tables
write.csv(x = round(del.AIC,2), file = 'INSIGHT_seasonal_AICs.csv')
write.csv(x = round(del.pdmAIC,2), file = 'INSIGHT_pandemic_AICs.csv')






