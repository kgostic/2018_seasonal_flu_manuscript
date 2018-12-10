rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')
load('H1ests_vv_OAS_2017-09-03.RData')
load('H3ests_vv_OAS2017-09-03.RData')

rm(wm.H3)
rm(wm.H1)
rm(wo.H1)
rm(wo.H3)
load('raw_weights.RData')



## Try to account for ascertainment
# Data


###### ------ No Imprinting ------ ######
# Simulate model best fit using OAS estimates
OAS.H3 = simulate.incidence(alpha.hat = H3.ests.vv.OAS$par['alpha.hat'], AA.hat = H3.ests.vv.OAS$par['AA.hat'], R0.hat = H3.ests.vv.OAS$par['R0.hat'], vv.hat = H3.ests.vv.OAS$par['vv.hat'], tau.hat = H3.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2015, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)

OAS.H1 = simulate.incidence(alpha.hat = H1.ests.vv.OAS$par['alpha.hat'], AA.hat = H1.ests.vv.OAS$par['AA.hat'], R0.hat = H1.ests.vv.OAS$par['R0.hat'], vv.hat = H1.ests.vv.OAS$par['vv.hat'], tau.hat = H1.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL)
# Plot using a blue line
use.years = rownames(H3.inputs)



###### ------ Set grids of age-specific and imprinting-specific modulators ------ ######
imp.grid = seq(1, 15, by = 1)
age.grid = 1:55
full.grid = expand.grid(age.grid, imp.grid)

age.grid = c(1, 5, 10, 20, 30, 40, seq(47, 55, by = 1))
HN.grid = expand.grid(age.grid, imp.grid, imp.grid)


# Define funciton to simulate and calculate the likelihood
ascertain_nll = function(age.par, HH, OAS.simulation, dat.in, wm){
  use.years = rownames(dat.in)
  dat.in = dat.in[use.years, ]
  wm = wm[use.years, ]
  # For valid years in which no cases are predicted, use the prediction from the previous valid year
  OAS.valid = OAS.simulation[use.years, ]*0 # Initialize
  for(yy in 1:length(use.years)){
    yr = use.years[yy]
    # If no cases occurred in the year of interest, go back in time until cases were predicted
    while(sum(OAS.simulation[yr, ]) <= 0){
      yr = as.character(as.numeric(yr) - 1)
    }
  OAS.valid[yy, ] = OAS.simulation[yr, ]    
  }
  wo = 1-wm
  # Assume risk increases linearly with age
  elderly.mat = matrix(c(rep(0, 65), 1:35*age.par/35+1), nrow = length(use.years), ncol = 100, byrow = T); colnames(elderly.mat) = 0:99
  adult.child.mat = elderly.mat*0; adult.child.mat[, as.character(0:64)] = 1
  observed =  OAS.valid * (wm + HH*wo) * (elderly.mat + adult.child.mat)
  
  #Normalize
  observed = observed/sum(observed)
  nll.yearly = vector('numeric')
  for(ii in 1:length(use.years)){
  nll.yearly[ii] = -dmultinom(x = as.numeric(dat.in[ii,]), prob = as.numeric(observed[ii, ]), log = TRUE)  }
  
 return( list(nll = sum(nll.yearly), observed = observed))
}



# Define funciton to simulate and calculate the likelihood for both H and N effects
ascertain_nll_HN = function(age.par, HH, NN, OAS.simulation, dat.in, wmm, wom, wmo, woo){
  use.years = rownames(dat.in)
  dat.in = dat.in[use.years, ]
  wmm = wmm[use.years, ]
  wmo = wmo[use.years, ]
  wom = wom[use.years, ]
  woo = woo[use.years, ]
  # For valid years in which no cases are predicted, use the prediction from the previous valid year
  OAS.valid = OAS.simulation[use.years, ]*0 # Initialize
  for(yy in 1:length(use.years)){
    yr = use.years[yy]
    # If no cases occurred in the year of interest, go back in time until cases were predicted
    while(sum(OAS.simulation[yr, ]) <= 0){
      yr = as.character(as.numeric(yr) - 1)
    }
    OAS.valid[yy, ] = OAS.simulation[yr, ]    
  }
  # Assume risk increases linearly with age
  elderly.mat = matrix(c(rep(0, 65), 1:35*age.par/35+1), nrow = length(use.years), ncol = 100, byrow = T); colnames(elderly.mat) = 0:99
  adult.child.mat = elderly.mat*0; adult.child.mat[, as.character(0:64)] = 1
  observed =  OAS.valid * (wmm + wmo*NN + wom*HH + woo*NN*HH) * (elderly.mat + adult.child.mat)
  
  #Normalize
  observed = observed/sum(observed)
  nll.yearly = vector('numeric')
  for(ii in 1:length(use.years)){
    nll.yearly[ii] = -dmultinom(x = as.numeric(dat.in[ii,]), prob = as.numeric(observed[ii, ]), log = TRUE)  }
  
  return( list(nll = sum(nll.yearly), observed = observed))
}



###################################################################
####                     H3N2
###################################################################

###### ------ Grid search for an AH, AN and ANH model ------ ######
# results.AH.H3 = mapply(FUN = ascertain_nll, age.par = full.grid[,1], HH = full.grid[,2], MoreArgs = list(OAS.simulation = OAS.H3, dat.in = H3.inputs, wm = w3.formatted), SIMPLIFY = TRUE)
# results.AN.H3 = mapply(FUN = ascertain_nll, age.par = full.grid[,1], HH = full.grid[,2], MoreArgs = list(OAS.simulation = OAS.H3, dat.in = H3.inputs, wm = w3.formatted+w2.formatted), SIMPLIFY = TRUE)
# results.AHN.H3 = mapply(FUN = ascertain_nll_HN, age.par = HN.grid[,1], HH = HN.grid[,2], NN = HN.grid[,3], MoreArgs = list(OAS.simulation = OAS.H3, dat.in = H3.inputs, wmm = w3.formatted, wom = w2.formatted, wmo = w2.formatted*0, woo = w1.formatted+wn.formatted), SIMPLIFY = TRUE)

load('Grid_results.RData')

###### ------ Grid search for an AH, AN and ANH model ------ ######
best.pars = function(res, grid = full.grid){
nll.vals = unlist(res[1, ])
sims = res[2, ]
which.min(nll.vals)
grid[which.min(nll.vals), ]
}

###### ------ Get the best parameters from each model ------ ######
best.pars(results.AH.H3)
best.pars(results.AN.H3)
best.pars(results.AHN.H3, grid = HN.grid)

valid = which(full.grid[,1] == 1)  ## Subset of no age effects
full.grid[valid, ]
which.min(results.AN.H3[1, valid]) ## Best imprinting only par

valid = which(full.grid[,2] == 1)  ## Subset of no imp effects
full.grid[valid, ]
which.min(results.AN.H3[1, valid]) ## Best age only par


###################################################################
####                     H1N1
###################################################################

###### ------ Grid search for an AH, AN and ANH model ------ ######
# results.AH.H1 = mapply(FUN = ascertain_nll, age.par = full.grid[,1], HH = full.grid[,2], MoreArgs = list(OAS.simulation = OAS.H1, dat.in = H1.inputs, wm = w1.formatted+w2.formatted), SIMPLIFY = TRUE)
# results.AN.H1 = mapply(FUN = ascertain_nll, age.par = full.grid[,1], HH = full.grid[,2], MoreArgs = list(OAS.simulation = OAS.H1, dat.in = H1.inputs, wm = w1.formatted), SIMPLIFY = TRUE)
# results.AHN.H1 = mapply(FUN = ascertain_nll_HN, age.par = HN.grid[,1], HH = HN.grid[,2], NN = HN.grid[,3], MoreArgs = list(OAS.simulation = OAS.H1, dat.in = H1.inputs, wmm = w1.formatted, wom = w2.formatted*0, wmo = w2.formatted, woo = w3.formatted+wn.formatted), SIMPLIFY = TRUE)


#save(results.AH.H1, results.AH.H3, results.AN.H1, results.AN.H3, results.AHN.H1, results.AHN.H3, file = 'Grid_results.RData')




###### ------ Grid search for an AH, AN and ANH model ------ ######
best.pars = function(res, grid = full.grid){
  nll.vals = unlist(res[1, ])
  sims = res[2, ]
  which.min(nll.vals)
  grid[which.min(nll.vals), ]
}

###### ------ Get the best parameters from each model ------ ######
best.pars(results.AH.H1)
best.pars(results.AN.H1)
best.pars(results.AHN.H1, grid = HN.grid)












#Plot data
cols = 'red'
plot(0:99, H3.inputs[1, ]/sum(H3.inputs[1, ]), col = cols, main = '', ylab = 'Fraction of cases', xlab = 'Case age', ylim = c(0, .09))
for(ii in 1:nrow(H3.inputs)){ points(0:99, H3.inputs[ii, ]/sum(H3.inputs[ii, ]), col = cols)}
mtext('H3N2', line =.5, font = 2)
points(0:99, colSums(H3.inputs)/sum(H3.inputs), pch = 16, col = 'pink2')
# Plot aggregate OAS
lines(0:99, colSums(OAS.H3)/sum(OAS.H3), col = 'white', lwd = 4)
# Age only
age.only = ascertain_nll(age.par = 55, HH = 1, OAS.simulation = OAS.H3, dat.in = H3.inputs, wm = w3.formatted)$observed
lines(0:99, colSums(age.only)/sum(age.only), col = 'yellow', lwd = 4)
legend(15, .085, c('Data - single season', 'Data - all seasons', 'Model - no imprinting', 'Model - age bias'), pch = c(1, 16, NA, NA), col = c('red', 'pink2', 'white', 'yellow'), bty = 'n', lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3))

 # Best 
bp = best.pars(results.AH.H3, grid = full.grid)
best.sim = results.AH.H3[2, ][[ which.min(results.AH.H3[1, ]) ]]
lines(0:99, colSums(best.sim)/sum(best.sim), col = 'purple', lwd = 4)
legend(15, .085, c('Data - single season', 'Data - all seasons', 'Model - no imprinting', 'Model - age bias', 'Model - age * imprinting bias'), pch = c(1, 16, NA, NA, NA), col = c('red', 'pink2', 'white', 'yellow', 'purple'), bty = 'n', lty = c(NA, NA, 1, 1, 1), lwd = c(NA, NA, 3, 3, 3))




## Plot grid space 
n.pars = rep(2, nrow(full.grid))
n.pars[which(full.grid[,1] == 1)] = n.pars[which(full.grid[,1] == 1)] - 1
n.pars[which(full.grid[,2] == 1)] = n.pars[which(full.grid[,2] == 1)] - 1

AIC = 2*n.pars + 2*nll.vals
which.min(AIC)
age.vals = matrix(full.grid[, 1], nrow = length(unique(full.grid[,1])))
imp.vals = matrix(full.grid[, 2], nrow = length(unique(full.grid[,1])))
del.AIC = AIC-min(AIC)
daic = matrix(del.AIC, nrow = length(unique(full.grid[,1])))

image(age.vals[,1], imp.vals[1,], daic, col = rainbow(40))
contour(age.vals[,1], imp.vals[1,], daic, add = TRUE, levels = c(seq(0, 100, by = 2), seq(200, 1800, by = 200)), col = 'black')











## H1N1



age.grid = seq(1, 10, by = 1)
imp.grid = seq(1, 4, by = 1)
full.grid = expand.grid(age.grid, imp.grid)


results.H = mapply(FUN = ascertain_nll, age.par = full.grid[,1], HH = full.grid[,2], MoreArgs = list(OAS.simulation = OAS.H1, dat.in = H1.inputs, wm = wm.H1), SIMPLIFY = TRUE)

nll.vals = unlist(results.H[1, ])
sims = results.H[2, ]

which.min(nll.vals)
full.grid[which.min(nll.vals), ]


cols = gray.colors(10)[-1]
plot(0:99, H1.inputs[1, ]/sum(H1.inputs[1, ]), ylim = c(0, .09), main = 'Imprinting', col = cols[1])
for(ii in 1:nrow(H1.inputs)){ points(0:99, H1.inputs[ii, ]/sum(H1.inputs[ii, ]), col = cols[ii])}
points(0:99, colSums(H1.inputs)/sum(H1.inputs), pch = 16)
lines(0:99, colSums(OAS.H1)/sum(OAS.H1), lwd = 2)
lines(0:99, colSums(sims[[1]])/sum(sims[[1]]), col = 'blue', lwd = 2)

n.pars = rep(2, nrow(full.grid))
n.pars[which(full.grid[,1] == 1)] = n.pars[which(full.grid[,1] == 1)] - 1
n.pars[which(full.grid[,2] == 1)] = n.pars[which(full.grid[,2] == 1)] - 1

AIC = 2*n.pars + 2*nll.vals
which.min(AIC)
age.vals = matrix(full.grid[, 1], nrow = length(unique(full.grid[,1])))
imp.vals = matrix(full.grid[, 2], nrow = length(unique(full.grid[,1])))
del.AIC = AIC-min(AIC)
daic = matrix(del.AIC, nrow = length(unique(full.grid[,1])))


image(age.vals[,1], imp.vals[1,], daic, col = rainbow(40))
contour(age.vals[,1], imp.vals[1,], daic, add = TRUE, levels = c(seq(0, 100, by = 2), seq(200, 1800, by = 200)), col = 'black')















## data only
par(bg = 'black', fg = 'white', col.main = 'white', col.axis = 'white', col.lab = 'white')
cols = 'red'
plot(0:99, H3.inputs[1, ]/sum(H3.inputs[1, ]), ylim = c(0, .09), main = 'H3N2', col = cols[1], ylab = 'Fraction of cases', xlab = 'Case age')
for(ii in 1:nrow(H3.inputs)){ points(0:99, H3.inputs[ii, ]/sum(H3.inputs[ii, ]), col = cols)}
points(0:99, colSums(H3.inputs)/sum(H3.inputs), pch = 16, col = 'pink3')
mtext(paste(sum(H3.inputs), ' cases'))
legend(45, .085, c('Single season', 'All seasons'), pch = c(1, 16), col = c('red', 'pink3'), bty = 'n')


cols = 'dodgerblue'
plot(0:99, H1.inputs[1, ]/sum(H1.inputs[1, ]), col = cols, main = 'H1N1', ylab = 'Fraction of cases', xlab = 'Case age', ylim = c(0, .09))
for(ii in 1:nrow(H1.inputs)){ points(0:99, H1.inputs[ii, ]/sum(H1.inputs[ii, ]), col = cols)}
points(0:99, colSums(H1.inputs)/sum(H1.inputs), pch = 16, col = 'skyblue')
mtext(paste(sum(H1.inputs), ' cases'))
legend(45, .085, c('Single season', 'All seasons'), pch = c(1, 16), col = c('dodgerblue', 'skyblue'), bty = 'n')



plot(0:99, colSums(H3.inputs)/sum(H3.inputs), pch = 16, col = 'red', ylim = c(0, .09), main = 'Observed data', ylab = 'Fraction of cases', xlab = 'Case age')
points(0:99, colSums(H1.inputs)/sum(H1.inputs), pch = 16, col = 'dodgerblue')
legend(45, .085, c('H1N1', 'H3N2'), col = c('dodgerblue', 'red'), pch = 16, bty = 'n')
mtext('All seasons')

