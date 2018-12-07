## Grid search

source('Baseline.R')
# 
pdf('Ht_visualization.pdf')

## Plot different outcomes
load('H1ests_vv_OAS_2017-09-01.RData')

# Simulate outcomes with just OAS pars
null = simulate.incidence(alpha.hat = H1.ests.vv.OAS$par['alpha.hat'], AA.hat = H1.ests.vv.OAS$par['AA.hat'], R0.hat = H1.ests.vv.OAS$par['R0.hat'], vv.hat = H1.ests.vv.OAS$par['vv.hat'], tau.hat = H1.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H1, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = FALSE)



# Simulate across a grid of H values
grid.H1 = lapply(seq(.2, 1, by = .1), FUN = function(Hm.val) simulate.incidence(alpha.hat = H1.ests.vv.OAS$par['alpha.hat'], AA.hat = H1.ests.vv.OAS$par['AA.hat'], R0.hat = H1.ests.vv.OAS$par['R0.hat'], vv.hat = H1.ests.vv.OAS$par['vv.hat'], tau.hat = H1.ests.vv.OAS$par['tau.hat'], Hm = Hm.val, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H1, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = TRUE))


# Plot grid
dev.off()
par(mfrow = c(1,1))
plot(0:99, colSums(H1.inputs)/sum(H1.inputs), ylim = c(0, .07), xlab = 'age', ylab = 'frac cases', main = 'H1N1', pch = 16)
cols = gray.colors(10)[-1]
for(ii in 1:nrow(H1.inputs)){ points(0:99, H1.inputs[ii, ]/sum(H1.inputs[ii, ]), col = cols[ii],cex = .7)}
use.rows = rownames(H1.inputs)
lines(0:99, colSums(null[use.rows, ])/sum(null[use.rows, ]))
cols = rainbow(15)[-c(1, 2, 4)]
for(ii in 1:length(grid.H1)){
  yy = grid.H1[[ii]][use.rows, ]
  lines(0:99, colSums(yy)/sum(yy), col = cols[ii])
}
legend('topright', c('data', 'prediction'), pch = c('o', '-'))
legend(50, .06, paste('Hm=', seq(.2, 1, by = .1)), lty = 1, col = cols, cex = .8)

#scaled.wm.H1 = wm.H1[rownames(H1.inputs), ]*rowSums(H1.inputs)/sum(H1.inputs)
#barplot(scaled.wm.H1, xlab = 'age', ylab = 'Group 1 Imprinting protection', main = 'G1 imprinting prob, scaled by case counts each year')


# Plot grid residuals
layout(matrix(c(1, 1, 1, 2), nrow = 4, ncol = 1))
obs = colSums(H1.inputs)/sum(H1.inputs)
cols = rainbow(15)
use.rows = rownames(H1.inputs)
yy = grid.H1[[1]][use.rows, ]
plot(0:99, obs-(colSums(yy)/sum(yy)), xlab = 'age', ylab = 'obs - exp', main = 'H1N1', col = cols[1], type = 'l')
abline(h = 0)
sum( (obs-(colSums(yy)/sum(yy)))^2)

for(ii in 2:length(grid.H1)){
  yy = grid.H1[[ii]][use.rows, ]
  lines(0:99, obs-colSums(yy)/sum(yy), col = cols[ii])
}
sum( (obs-(colSums(yy)/sum(yy)))^2)
legend('topright', c('data', 'prediction'), pch = c('o', '-'))
legend(50, -.01, paste('Hm=', seq(0, 1, by = .1)), lty = 1, col = cols, cex = .8)

scaled.wm.H1 = wm.H1[rownames(H1.inputs), ]*rowSums(H1.inputs)/sum(H1.inputs)
barplot(scaled.wm.H1, xlab = 'age', ylab = 'Group 1 Imprinting protection', main = 'G1 imprinting prob, scaled by case counts each year')





load('H3ests_vv_OAS2017-09-02.RData')
null = simulate.incidence(alpha.hat = H3.ests.vv.OAS$par['alpha.hat'], AA.hat = H3.ests.vv.OAS$par['AA.hat'], R0.hat = H3.ests.vv.OAS$par['R0.hat'], vv.hat = H3.ests.vv.OAS$par['vv.hat'], tau.hat = H3.ests.vv.OAS$par['tau.hat'], Hm = 1, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H1, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = FALSE)



# Simulate across a grid of H values
grid.H3 = lapply(seq(.2, 1, by = .1), FUN = function(Hm.val) simulate.incidence(alpha.hat = H3.ests.vv.OAS$par['alpha.hat'], AA.hat = H3.ests.vv.OAS$par['AA.hat'], R0.hat = H3.ests.vv.OAS$par['R0.hat'], vv.hat = H3.ests.vv.OAS$par['vv.hat'], tau.hat = H3.ests.vv.OAS$par['tau.hat'], Hm = Hm.val+((1-Hm.val)*.2), demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H3, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = TRUE))


# Plot grid
dev.off()
par(mfrow = c(1,1))
plot(0:99, colSums(H3.inputs)/sum(H3.inputs), xlab = 'age', ylab = 'frac cases', main = 'H3N2', ylim = c(0, .06), pch = 16)
cols = gray.colors(10)
for(ii in 1:nrow(H3.inputs)){ points(0:99, H3.inputs[ii, ]/sum(H3.inputs[ii, ]), col = cols[ii], cex = .7)}
use.rows = rownames(H3.inputs)
lines(0:99, colSums(null[use.rows, ])/sum(null[use.rows, ]))
cols = rainbow(15)[-c(1, 2, 4)]
for(ii in length(grid.H3):1){
  yy = grid.H3[[ii]][use.rows, ]
  lines(0:99, colSums(yy)/sum(yy), col = cols[ii])
  #Sys.sleep(3)
  #print(ii)
}
legend('topright', c('data', 'prediction'), pch = c('o', '-'))
legend(50, .06, paste('Hm=', seq(.2, 1, by = .1)), lty = 1, col = cols, cex = .8)

scaled.wm.H3 = wm.H3[rownames(H3.inputs), ]*rowSums(H3.inputs)/sum(H3.inputs)
barplot(scaled.wm.H3, xlab = 'age', ylab = 'Group 2 Imprinting protection', main = 'G2 imprinting prob, scaled by case counts each year')






load('H1ests_vv_OAS_H_2017-08-16.RData')
load('H3ests_vv_OAS_H_2017-08-16.RData')
# Plot grid residuals
layout(matrix(c(1, 1, 1, 2), nrow = 4, ncol = 1))
obs = colSums(H3.inputs)/sum(H3.inputs)
cols = rainbow(15)
use.rows = rownames(H3.inputs)
yy = grid.H3[[1]][use.rows, ]
plot(0:99, obs-(colSums(yy)/sum(yy)), xlab = 'age', ylab = 'obs - exp', main = 'H3N2', col = cols[1], type = 'l')
abline(h = 0)
sum( (obs-(colSums(yy)/sum(yy)))^2)

for(ii in 2:length(grid.H3)){
  yy = grid.H3[[ii]][use.rows, ]
  lines(0:99, obs-colSums(yy)/sum(yy), col = cols[ii])
}
sum( (obs-(colSums(yy)/sum(yy)))^2)
legend('topright', c('data', 'prediction'), pch = c('o', '-'))
legend(50, -.01, paste('Hm=', seq(0, 1, by = .1)), lty = 1, col = cols, cex = .8)

scaled.wm.H3 = wm.H3[rownames(H3.inputs), ]*rowSums(H3.inputs)/sum(H3.inputs)
barplot(scaled.wm.H3, xlab = 'age', ylab = 'Group 2 Imprinting protection', main = 'G2 imprinting prob, scaled by case counts each year')

dev.off()






#### Do a grid search!!
# library(parallel)
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
# ## export relevant variables
# clusterExport(cl, c('epi_final_size', 'simulate.incidence', 'demog', 'H1.clusters', 'H1.inputs', 'fraction.vaccinated', 'wm.H1', 'wo.H1', 'waifw.disaggregate', 'get_first_exposure_probs', 'get_last_exposure_probs', 'H1.ests.vv.OAS', 'nll'))
# 
# grid.vals.H1 = parSapply(cl, seq(0, .975, by = .025), FUN = function(Hm.val) nll(pars = c(H1.ests.vv.OAS$par, Hm = Hm.val), demog = demog, clusters = H1.clusters, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H1, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, dat.in = H1.inputs))
# 
# stopCluster(cl)
# 
# 
# 
# 
# #### Do a grid search!!
# library(parallel)
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores)
# ## export relevant variables
# clusterExport(cl, c('epi_final_size', 'simulate.incidence', 'demog', 'H3.clusters', 'H3.inputs', 'fraction.vaccinated', 'wm.H3', 'wo.H3', 'waifw.disaggregate', 'get_first_exposure_probs', 'get_last_exposure_probs', 'H3.ests.vv.OAS', 'nll'))
# 
# grid.vals.H3 = parSapply(cl, seq(0, .975, by = .001), FUN = function(Hm.val) nll(pars = c(H3.ests.vv.OAS$par, Hm = Hm.val), demog = demog, clusters = H3.clusters, st.yr = 1968, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, dat.in = H3.inputs))
# 
# stopCluster(cl)
# 
# save(grid.vals.H1, grid.vals.H3, file = 'Grid_search_Ht_vals.RData')

load('Grid_search_Ht_vals.RData')


pdf('Slice_nlls.pdf')
plot(seq(0, .975, by = .025), grid.vals.H1[1,], ylab = 'neg log lk', xlab = 'Hm value', main = 'H1N1 fits', type = 'b', pch = as.character(grid.vals.H1[2,]))
plot(seq(0, .975, by = .001), grid.vals.H3[1,], ylab = 'neg log lk', xlab = 'Hm value', main = 'H3N2 fits', type = 'b', pch = as.character(grid.vals.H3[2,]))
dev.off()


