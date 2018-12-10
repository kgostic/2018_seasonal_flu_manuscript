## Plot model comparison results
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')

# Load H1 model fits
load('H1ests_2017-08-01.RData'); H1.ests # baseline
load('H1ests_vv2017-08-01.RData'); H1.ests.vv # baseline + vaccination
load('H1ests_vv_OAS2017-08-01.RData'); H1.ests.vv.OAS # baseline + vaccination + OAS
load('H1ests_vv_OAS_H2017-08-01.RData'); H1.ests.vv.OAS.H # baseline + vaccination + OAS + HA imprinting effects on case ascertainment
load('H1ests_vv_OAS_Ht2017-08-02.RData'); H1.ests.vv.OAS.Ht # baseline + vaccination + OAS + HA imprinting effects on transmission
load('H1ests_vv_OAS_Hs2017-08-02.RData'); H1.ests.vv.OAS.Hs # baseline + vaccination + OAS + HA imprinting effects on transmission subtype-specific
load('H1ests_vv_LAS2017-08-01.RData'); H1.ests.vv.LAS # baseline + vaccination + LAS

# Load H3 model fits
load('H3ests_2017-07-29.RData'); H3.ests # baseline
load('H3ests_vv2017-08-01.RData'); H3.ests.vv # baseline + vaccination
load('H3ests_vv_OAS2017-08-02.RData'); H3.ests.vv.OAS # baseline + vaccination + OAS
load('H3ests_vv_OAS_H2017-08-05.RData'); H3.ests.vv.OAS.H # baseline + vaccination + OAS + HA imprinting
load('H3ests_vv_OAS_Ht2017-08-03.RData'); H3.ests.vv.OAS.Ht # baseline + vaccination+ OAS + HA imprinting
load('H3ests_vv_LAS2017-08-02.RData'); H3.ests.vv.LAS # baseline + vaccination + LAS


H1.results.table = H3.results.table = matrix(NA, nrow = 4, ncol = 6, dimnames = list(c('base', 'vax', 'vax_OAS', 'vax_OAS_HA'), names(H1.ests.vv.OAS.Ht$par)))


# Fill in results table
H1.results.table[1, names(H1.ests$par)] = H1.ests$par
H1.results.table[2, names(H1.ests.vv$par)] = H1.ests.vv$par
H1.results.table[3, names(H1.ests.vv.OAS$par)] = H1.ests.vv.OAS$par
H1.results.table[4, names(H1.ests.vv.OAS.Ht$par)] = H1.ests.vv.OAS.Ht$par

H3.results.table[1, names(H3.ests$par)] = H3.ests$par
H3.results.table[2, names(H3.ests.vv$par)] = H3.ests.vv$par
H3.results.table[3, names(H3.ests.vv.OAS$par)] = H3.ests.vv.OAS$par
H3.results.table[4, names(H3.ests.vv.OAS.Ht$par)] = H3.ests.vv.OAS.Ht$par

write.csv(x = H1.results.table, file = 'H1_par_ests.csv')
write.csv(x = H3.results.table, file = 'H3_par_ests.csv')

# Load  model predictions
load('Baseline_predictions.RData')
load('Baseline_vaccination_predictions.RData')
load('Baseline_vaccination_OAS_predictions.RData')
load('Baseline_vaccination_OAS_H_predictions.RData')
load('Baseline_vaccination_OAS_Ht_predictions.RData')
#load('Baseline_vaccination_LAS_predictions.RData')


# Import demographic data from the USA for all available years (1980-2020)
demog.raw = read.csv('USA_demography_1980_2020_formatted.csv', header = T)
# Drop individuals over age 99
demog = demog.raw[,-101]/rowSums(demog.raw[,-101]) # Normalize
rownames(demog) = 1980:2020; colnames(demog) = 0:99


## Calculate AIC
nlls.H1 = c(base = H1.ests$value, 
          V = H1.ests.vv$value, 
          V.OAS = H1.ests.vv.OAS$value, 
          V.LAS = H1.ests.vv.LAS$value, 
          V.OAS.H = H1.ests.vv.OAS.H$value,
          V.OAS.Ht = H1.ests.vv.OAS.Ht$value,
          V.OAS.Hs = H1.ests.vv.OAS.Hs$value)
n.pars = c(3, 4, 5, 5, 6, 6, 6)
AIC.H1 = 2*n.pars+2*nlls.H1; AIC.H1
del.AIC.H1 = sort(AIC.H1 - min(AIC.H1))

nlls.H3 = c(base = H3.ests$value, 
          V = H3.ests.vv$value,
          V.OAS = H3.ests.vv.OAS$value, 
          V.LAS = H3.ests.vv.LAS$value,
          V.OAS.H = H3.ests.vv.OAS.H$value,
          V.OAS.Ht = H3.ests.vv.OAS.Ht$value)
n.pars = c(3, 4, 5, 5, 6, 6)
AIC.H3 = 2*n.pars+2*nlls.H3; AIC.H3
del.AIC.H3 = sort(AIC.H3 - min(AIC.H3))

del.AIC.H1
del.AIC.H3

############################################################
#____________________ Load NCBI data ______________________#
#
############################################################
H1_NCBI = read.csv('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/NCBI_H1_Data.csv', stringsAsFactors = FALSE)
H3_NCBI = read.csv('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/NCBI_H3_Data.csv', stringsAsFactors = FALSE)

## Load Arizona data
H1_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H1.csv')
H3_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H3.csv')


################# 
reformat.data = function(dat.in){
  yrs = unique(dat.in$Year)
  ages = 0:99
  out.mat = matrix(NA, nrow = length(yrs), ncol = length(ages), dimnames = list(yrs, ages))
  for(yy in 1:length(yrs)){
    for(aa in 1:length(ages)){
      out.mat[yy, aa] = sum(dat.in$Year == yrs[yy] & dat.in$Age == ages[aa])
    }
  }
  out.mat
}

## Reformatted data inputs
H1.inputs = reformat.data(H1_AZ)
H3.inputs = reformat.data(H3_AZ)

H1.inputs = H1.inputs[-which(rowSums(H1.inputs) < 100), ]
H3.inputs = H3.inputs[-which(rowSums(H3.inputs) < 100), ]





pdf('All_model_fits.pdf')
### PLOT
dat.years = rownames(H1.inputs)
par(mfrow = c(4, 2), mar = c(3, 2, 2, 2))
for(ii in 1:nrow(H1.inputs)){
  yr = rownames(H1.inputs)[ii]
  plot(0:99, H1.inputs[ii, ]/rowSums(H1.inputs)[ii], xlab = 'Age', ylab = 'Fraction of cases', main = paste('H1N1', rownames(H1.inputs)[ii]), pch = 16)
  points(0:99, H1.sim.baseline.vv[yr, ], col = 'blue', cex = .5)
  points(0:99, H1.sim.baseline.vv.OAS[yr, ], col = 'green3', cex = .5)
  points(0:99, H1.sim.baseline.vv.OAS.H[yr, ], col = 'violet', cex = .5)
  legend('topright', c('Data', 'Cross immunity + vaccination', '...+ OAS', '...+ OAS + HA Imprinting'), lty = c(16, 1, 1, 1), col = c('black', 'blue', 'green3', 'violet'), bty = 
           'n', cex = .7) 
}
plot(0:99, colSums(H1.inputs)/sum(H1.inputs), main = 'Overall', pch = 16)
points(0:99, colSums(H1.sim.baseline.vv)/sum(H1.sim.baseline.vv), col = 'blue', cex = .5)
points(0:99, colSums(H1.sim.baseline.vv.OAS)/sum(H1.sim.baseline.vv.OAS), col = 'green3', cex = .5)
points(0:99, colSums(H1.sim.baseline.vv.OAS.H)/sum(H1.sim.baseline.vv.OAS.H), col = 'violet', cex = .5)
legend('topright', c('Data', 'Cross immunity + vaccination', '...+ OAS', '...+ OAS + HA Imprinting'), pch = c(16, 1, 1, 1), col = c('black', 'blue', 'green3', 'violet'), bty = 
         'n', cex = .7) 

par(mfrow = c(1,1), mar = c(4,4,4,4))
plot(0:99, colSums(H1.inputs)/sum(H1.inputs), main = 'H1N1 Overall', pch = 16, ylab = 'fraction of cases in age group', xlab = 'age group')
points(0:99, colSums(H1.sim.baseline.vv)/sum(H1.sim.baseline.vv), col = 'blue', cex = .5)
points(0:99, colSums(H1.sim.baseline.vv.OAS)/sum(H1.sim.baseline.vv.OAS), col = 'green3', cex = .5)
points(0:99, colSums(H1.sim.baseline.vv.OAS.H)/sum(H1.sim.baseline.vv.OAS.H), col = 'violet', cex = .5)
legend('topright', c('Data', 'Cross immunity + vaccination', '...+ OAS', '...+ OAS + HA Imprinting'), pch = c(16, 1, 1, 1), col = c('black', 'blue', 'green3', 'violet'), bty = 
         'n', cex = 1) 


#H3N2 
dat.years = rownames(H3.inputs)
par(mfrow = c(4, 2), mar = c(3, 2, 2, 2))
yrs = H3.inputs[-c(4, 7, 10)]
for(ii in c(1,2,3,5,6,8,9)){
#for(ii in 1:10){
  yr = yrs[ii]
  plot(0:99, H3.inputs[ii, ]/rowSums(H3.inputs)[ii], xlab = 'Age', ylab = 'Fraction of cases', main = paste('H3N1', rownames(H3.inputs)[ii]), pch = 16)
  points(0:99, H3.sim.baseline.vv[yr, ], col = 'blue', cex = .7)
  points(0:99, H3.sim.baseline.vv.OAS[yr, ], col = 'green3', cex = .7)
  points(0:99, H3.sim.baseline.vv.OAS.H[yr, ], col = 'violet', cex = .7)
  legend('topright', c('Data', 'Cross immunity + vaccination', '...+ OAS', '...+ OAS + HA Imprinting'), pch = c(16, 1, 1, 1), col = c('black', 'blue', 'green3', 'violet'), bty = 
           'n', cex = .7) 
}
plot(0:99, colSums(H3.inputs)/sum(H3.inputs), main = 'Overall', pch = 16)
points(0:99, colSums(H3.sim.baseline.vv)/sum(H3.sim.baseline.vv), col = 'blue', cex = .7)
points(0:99, colSums(H3.sim.baseline.vv.OAS)/sum(H3.sim.baseline.vv.OAS), col = 'green3', cex = .7)
points(0:99, colSums(H3.sim.baseline.vv.OAS.H)/sum(H3.sim.baseline.vv.OAS.H), col = 'violet', cex = .7)
legend('topright', c('Data', 'Cross immunity + vaccination', '...+ OAS', '...+ OAS + HA Imprinting'), pch = c(16, 1, 1, 1), col = c('black', 'blue', 'green3', 'violet'), bty = 
         'n', cex = .7) 

par(mfrow = c(1,1))
plot(0:99, colSums(H3.inputs)/sum(H3.inputs), main = 'H3N2 Overall', pch = 16, ylab = 'fraction of cases in age group', xlab = 'age group')
points(0:99, colSums(H3.sim.baseline.vv)/sum(H3.sim.baseline.vv), col = 'blue', cex = .7)
points(0:99, colSums(H3.sim.baseline.vv.OAS)/sum(H3.sim.baseline.vv.OAS), col = 'green3', cex = .7)
points(0:99, colSums(H3.sim.baseline.vv.OAS.H)/sum(H3.sim.baseline.vv.OAS.H), col = 'violet', cex = .7)
legend('topright', c('Data', 'Cross immunity + vaccination', '...+ OAS', '...+ OAS + HA Imprinting'), pch = c(16, 1, 1, 1), col = c('black', 'blue', 'green3', 'violet'), bty = 
         'n', cex = 1) 
dev.off()



## Likelhood ratio tests
## Compute LR statistic
## See pages 191-192 of Bolker, 2008
load('H1ests_vv_OAS2017-08-01.RData')
load('H3ests_vv_OAS2017-08-02.RData')
load('Profs_H1_Ht2017-08-04.RData')
load('Profs_H3_Ht.RData')

## H1N1 fits. Test Ht vs. OAS only.
restricted.nll = H1.ests.vv.OAS$value
alternative.nll = H1.OAS.Ht.profs[[34]]$value
LR.stat = 2*(restricted.nll - alternative.nll)
## LR stat is distributed chi square with 1 df
pchisq(LR.stat, df = 1)
plot(seq(0, 5, by = .1), pchisq(seq(0, 5, by = .1), 1), xlab = 'LR stat', ylab = 'Cumulative density', bty = 'n')
abline(v = LR.stat, col = 'red')
abline(h = .95); text(3.6, .97, '.95')
abline(v = qchisq(.95, 1))



## H3N2 fits. Test Ht vs. OAS only.
restricted.nll = H3.ests.vv.OAS$value
alternative.nll = H3.OAS.Ht.profs[[16]]$value
LR.stat = 2*(restricted.nll - alternative.nll)
## LR stat is distributed chi square with 1 df
pchisq(LR.stat, df = 1)
plot(seq(0, 150, by = .1), pchisq(seq(0, 150, by = .1), 1), xlab = 'LR stat', ylab = 'Cumulative density', bty = 'n')
abline(v = LR.stat, col = 'red')
abline(h = .95); text(0, .97, '.95')
abline(v = qchisq(.95, 1))
