## alt model
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')
load('H1ests_vv_OAS_2017-08-16.RData')
load('H3ests_vv_OAS2017-08-16.RData')

H3.base = simulate.incidence(alpha.hat = H3.ests.vv.OAS$par['alpha.hat'], AA.hat = H3.ests.vv.OAS$par['AA.hat'], vv.hat = H3.ests.vv.OAS$par['vv.hat'], tau.hat = H3.ests.vv.OAS$par['tau.hat'], R0.hat = H3.ests.vv.OAS$par['R0.hat'], Hm = 1, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = TRUE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)

H1.base = simulate.incidence(alpha.hat = H1.ests.vv.OAS$par['alpha.hat'], AA.hat = H1.ests.vv.OAS$par['AA.hat'], vv.hat = H1.ests.vv.OAS$par['vv.hat'], tau.hat = H1.ests.vv.OAS$par['tau.hat'], R0.hat = H1.ests.vv.OAS$par['R0.hat'], Hm = 1, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = TRUE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)


wm = wm.H3[rownames(H3.inputs), ]
wo = 1-wm




## Define nll function

nll_alt = function(pars, dat.in, plot = FALSE, wm.in, p0, clust.in){

  aa = pars['aa']  
  if(aa < 0 | aa > 1) { stop('aa must be between 0 and 1')} 

    #___________________ calculate neg log lk ____________________ #  
    fit.years = rownames(dat.in)
    fit.clust = clust.in[as.character(clust.in$year) %in% fit.years, ]
    wm = wm.in[as.character(fit.years), ]
    wo = 1-wm
    wo = wo/rowSums(wo)
    p.hat = p0[fit.years, ]*aa + (1-aa)*wo # aa gives the weight placed on seasonal, 1-aa gives the weight placed on imprinting
    
    
    for(ii in 1:nrow(dat.in)){
    plot(0:99,dat.in[ii,]/sum(dat.in[ii,]), main = fit.clust$jump.year[ii])
    lines(0:99, p0[fit.years[ii], ], col = 'red', cex = .5)
    lines(0:99, p.hat[ii,], col = 'blue', cex = .5)
    }
    
    nll.yearly = vector('numeric', length(fit.years))
    counter = 0
    for(ii in 1:length(fit.years)){
      if( sum(p.hat[ii, ]) == 0 ){
        #print(ii)
        counter = counter + 1
        nll.yearly[ii] = penalty  # Penalize if no outbreak occurs in the year of interest
      }else{
        nll.yearly[ii] = -dmultinom(x = as.numeric(dat.in[ii,]), prob = as.numeric(p.hat[ii,]), log = TRUE)
      }
    }
  sum(nll.yearly)
  #nll.yearly
}




## Test full data set
H3.alt = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H3.inputs, wm.in = wm.H3, p0 = H3.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = H3.clusters); H3.alt
H1.alt = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H1.inputs, wm.in = wm.H1, p0 = H1.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = H1.clusters); H1.alt
load('Group1_weights.RData')
H1.subtype.alt = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H1.inputs, wm.in = w1.formatted, p0 = H1.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = H1.clusters); H1.subtype.alt

## Test H1 repeat vs. cluster jump years only
H1.dat.yrs = rownames(H1.inputs)
H1.dat.clust = H1.clusters[as.character(H1.clusters$year) %in% H1.dat.yrs, ]
H1.jump = H1.inputs['2009', ]
jump.clust.H1 = H1.clusters[33, ]
H1.no.jump = H1.inputs[c('2008', '2010', '2011', '2013'), ]
nojump.clust.H1 = H1.clusters[c(32, 34, 35, 37), ]

H1.alt.jump = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H1.jump, wm.in = wm.H1, p0 = H1.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = jump.clust.H1); H1.alt.jump
H1.alt.nojump = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H1.no.jump, wm.in = wm.H1, p0 = H1.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = nojump.clust.H1); H1.alt.nojump


## Test H3 repeat vs. cluster jump years only
H3.dat.yrs = rownames(H3.inputs)
H3.dat.clust = H3.clusters[as.character(H3.clusters$year) %in% H3.dat.yrs, ]
H3.jump = H3.inputs[c('2005', '2011'), ]
jump.clust.H3 = H3.clusters[c(38, 44), ]
H3.no.jump = H3.inputs[c('2008', '2010', '2012', '2014'), ]
nojump.clust.H3 = H3.clusters[c(41, 43, 45, 47), ]
H3.alt.jump = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H3.jump, wm.in = wm.H3, p0 = H3.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = jump.clust.H3); H3.alt.jump
H3.alt.nojump = optim(par = c('aa' = .8), fn = nll_alt, dat.in = H3.no.jump, wm.in = wm.H3, p0 = H3.base, method = 'L-BFGS-B', lower = 0, upper = 1, clust.in = nojump.clust.H3); H3.alt.nojump



## No real difference between jump and non-jump years!





# Plot data
par(mfrow = c(3,2), mar = c(4, 4, 2.5, 1), bg = 'black', fg = 'white', col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white')
dat.cols = colorRampPalette(c('white', 'black'))(10)


# Plot predictions
use.ests = H3.base[rownames(H3.inputs), ]
use.ests = use.ests/rowSums(use.ests)
cols1 = colorRampPalette(c('red', 'orange'))(6)

a = H3.alt$par
use.ests.5 = (use.ests*a) + (wo/rowSums(wo)*(1-a))
use.ests.5 = use.ests.5/rowSums(use.ests.5)
cols2 = colorRampPalette(c('blue', 'cyan'))(6)


plot(0:99, H3.inputs[1, ]/sum(H3.inputs[1,]), ylim = c(0, .07), main = rownames(H3.inputs)[1], ylab = 'Case frequency', xlab = 'age')
lines(0:99, use.ests[ii,], col = cols1[1])
lines(0:99, use.ests.5[ii,], col = cols2[1], lty = 1)
legend('topright', c('data', 'Age effects', 'Age + Birth Year effects'), pch = c(1, NA, NA), lty = c(NA, 1, 1), col = c('black', 'orangered', 'dodgerblue'), cex = .7, bty = 'n')

for(ii in 2:nrow(H3.inputs)){
  plot(0:99, H3.inputs[ii,]/sum(H3.inputs[ii,]), col = dat.cols[ii], main = rownames(H3.inputs)[ii], ylab = 'Case frequency', xlab = 'age', ylim = c(0, .07))
  lines(0:99, use.ests[ii,], col = cols1[ii])
  lines(0:99, use.ests.5[ii,], col = cols2[ii], lty = 1)
}








# Plot data
use.ests = H1.base[rownames(H1.inputs), ]
use.ests = use.ests/rowSums(use.ests)
a = H1.alt$par
use.ests.5 = (use.ests*a) + (wo/rowSums(wo)*(1-a))
use.ests.5 = use.ests.5/rowSums(use.ests.5)

par(mfrow = c(3,2), mar = c(3, 3, 1, 1), bg = 'black', fg = 'white', col = 'white', col.axis = 'white', col.lab = 'white', col.main = 'white')
for(ii in 1:nrow(H1.inputs)){
  plot(0:99, H1.inputs[ii,]/sum(H1.inputs[ii,]), col = dat.cols[ii], main = rownames(H1.inputs)[ii], ylab = '', xlab = '', ylim = c(0, .065))
  lines(0:99, use.ests[ii,], col = cols1[ii], lwd = 2)
  lines(0:99, use.ests.5[ii,], col = cols2[ii], lty = 2, lwd = 2)
  if(ii %% 2 == 1){
    mtext(text = 'Case frequency', side = 2, line = 2, cex = .8)
  }
  if(ii %in% c(5, 4)){
    mtext(text = 'Age', side = 1, line = 2, cex = .8)
  }
}
legend('topright', c('data', 'Age effects', 'Age + Birth Year effects'), pch = c(1, NA, NA), lty = c(NA, 1, 2), col = c('black', 'orangered', 'dodgerblue'), cex = .8, bty = 'n')

# Plot predictions
use.ests = H1.base[rownames(H1.inputs), ]
use.ests = use.ests/rowSums(use.ests)
cols = colorRampPalette(c('red', 'orange'))(6)
for(ii in 1:nrow(H1.inputs)){
  lines(0:99, use.ests[ii,], col = cols[ii])
}

a = H1.alt$par
use.ests.5 = (use.ests*a) + (wo/rowSums(wo)*(1-a))
use.ests.5 = use.ests.5/rowSums(use.ests.5)
cols = colorRampPalette(c('blue', 'cyan'))(6)
for(ii in 1:nrow(H1.inputs)){
  lines(0:99, use.ests.5[ii,], col = cols[ii], lty = 2)
}

legend('topright', c('data', 'Age effects', 'Age + Birth Year effects'), pch = c(1, NA, NA), lty = c(NA, 1, 2), col = c('black', 'orangered', 'dodgerblue'))




### Are results different in years with cluster jumps?!
