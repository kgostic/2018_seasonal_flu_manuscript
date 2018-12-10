## Grid search

source('Baseline.R')

## Plot different outcomes
load('H1ests_vv_OAS_2017-09-02.RData')

# Simulate outcomes with just OAS pars
null = simulate.incidence(alpha.hat = H1.ests.vv.OAS$par['alpha.hat'], AA.hat = H1.ests.vv.OAS$par['AA.hat'], R0.hat = H1.ests.vv.OAS$par['R0.hat'], vv.hat = H1.ests.vv.OAS$par['vv.hat'], tau.hat = H1.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H1, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = FALSE)


load('H1ests_vv_OAS_H_2017-09-02.RData')
best = simulate.incidence(alpha.hat = H1.ests.vv.OAS.H$par['alpha.hat'], AA.hat = H1.ests.vv.OAS.H$par['AA.hat'], R0.hat = H1.ests.vv.OAS.H$par['R0.hat'], vv.hat = H1.ests.vv.OAS.H$par['vv.hat'], tau.hat = H1.ests.vv.OAS.H$par['tau.hat'], Hm = H1.ests.vv.OAS.H$par['Hm'], demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H1, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = FALSE)


load('H1_Profs_coarse2017-09-03.RData')

H1.sims = vector('list', length(H1.OAS.Ht.profs))
Hms = seq(.4, 1, by = .1)
for(ii in 3:length(H1.OAS.Ht.profs)){
pp = H1.OAS.Ht.profs[[ii]]$par
H1.sims[[ii]] = simulate.incidence(alpha.hat = pp['alpha.hat'], AA.hat = pp['AA.hat'], R0.hat = pp['R0.hat'], tau.hat = pp['tau.hat'], vv.hat = pp['vv.hat'], Hm = Hms[ii], demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H1, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = FALSE)
print(ii)
print(Hms[ii])
}
Hms = Hms[-c(1,2)]
H1.sims = H1.sims[-c(1,2)]



dev.off()

# Plot grid
cols = rainbow(10)[-c(1, 3,4,6)]
layout(matrix(c(1, 1, 1, 2), nrow = 4, ncol = 1))
plot(0:99, colSums(H1.inputs)/sum(H1.inputs), ylim = c(0, .07), xlab = 'age', ylab = 'frac cases', main = 'H1N1')
use.rows = rownames(H1.inputs)
lines(0:99, colSums(null[use.rows, ])/sum(null[use.rows, ]), col = 'black')
#lines(0:99, colSums(best[use.rows, ])/sum(best[use.rows, ]), col = 'blue')

for(ii in 5:1){
lines(0:99, colSums(H1.sims[[ii]][use.rows, ])/sum(H1.sims[[ii]][use.rows, ]), col = cols[ii])
}


legend('topright', legend = paste('Hm =', Hms), lty = 1, col = cols[1:length(Hms)])








load('H3_Profs_coarse2017-09-03.RData')
Hms = seq(.4, 1, by = .1)
H3.sims = vector('list', length(H3.OAS.Ht.profs))
for(ii in 1:length(H3.OAS.Ht.profs)){
  pp = H3.OAS.Ht.profs[[ii]]$par
  H3.sims[[ii]] = simulate.incidence(alpha.hat = pp['alpha.hat'], AA.hat = pp['AA.hat'], R0.hat = pp['R0.hat'], tau.hat = pp['tau.hat'], vv.hat = pp['vv.hat'], Hm = Hms[ii], demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, wm.in = wm.H3, plot = FALSE, vaccinate = TRUE, WAIFW = waifw.disaggregate, OAS = TRUE, imprinting = FALSE)
}



# Plot grid
cols = rainbow(10)[-c(1, 3,4,6)]
layout(matrix(c(1, 1, 1, 2), nrow = 4, ncol = 1))
plot(0:99, colSums(H3.inputs)/sum(H3.inputs), ylim = c(0, .07), xlab = 'age', ylab = 'frac cases', main = 'H3N1')
use.rows = rownames(H3.inputs)
for(ii in 1:length(H3.sims)){
  lines(0:99, colSums(H3.sims[[ii]][use.rows, ])/sum(H3.sims[[ii]][use.rows, ]), col = cols[ii])
}
legend('topright', legend = paste('Hm =', Hms), lty = 1, col = cols[1:length(Hms)])

