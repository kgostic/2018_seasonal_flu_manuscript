setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')

## Compare model fits with and without imprinting
load('H1ests_vv_OAS_H_2017-09-04.RData')
sim.H = simulate.incidence(alpha.hat = H1.ests.vv.OAS.H$par['alpha.hat'], AA.hat = H1.ests.vv.OAS.H$par['AA.hat'], R0.hat = H1.ests.vv.OAS.H$par['R0.hat'], vv.hat = H1.ests.vv.OAS.H$par['vv.hat'], tau.hat = H1.ests.vv.OAS.H$par['tau.hat'], Hm = H1.ests.vv.OAS.H$par['Hm'], demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE,imprinting = TRUE, wm.in = wm.H1)

sim.adams.pars = simulate.incidence(alpha.hat = 0.25, AA.hat = 0.34, R0.hat = 1.24, vv.hat =0.06, tau.hat = 0.93, Hm = NULL, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H1)

load('H1ests_vv_OAS_2017-09-03.RData')
sim.OAS = simulate.incidence(alpha.hat = H1.ests.vv.OAS$par['alpha.hat'], AA.hat = H1.ests.vv.OAS$par['AA.hat'], R0.hat = H1.ests.vv.OAS$par['R0.hat'], vv.hat = H1.ests.vv.OAS$par['vv.hat'], tau.hat = H1.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H1)


# Plot data
cols = 'dodgerblue'
plot(0:99, H1.inputs[1, ]/sum(H1.inputs[1, ]), col = cols, main = '', ylab = 'Fraction of cases', xlab = 'Case age', ylim = c(0, .07))
for(ii in 1:nrow(H1.inputs)){ points(0:99, H1.inputs[ii, ]/sum(H1.inputs[ii, ]), col = cols)}
mtext('H1N1', line =.5, font = 2)
points(0:99, colSums(H1.inputs)/sum(H1.inputs), pch = 16, col = 'skyblue')
# Plot aggregate OAS
lines(0:99, colSums(sim.OAS)/sum(sim.OAS), col = 'white', lwd = 4)
lines(0:99, colSums(sim.OAS)/sum(sim.OAS), col = 'white', lwd = 2.5)
#lines(0:99, colSums(sim.adams.pars)/sum(sim.adams.pars), col = 'red')
legend(15, .07, c('Data - single season', 'Data - all seasons', 'Model - no imprinting'), pch = c(1, 16, NA), col = c('dodgerblue', 'skyblue', 'white'), bty = 'n', lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3))
#legend(40, .07, 'No imprinting', bty = 'n', lty = 1, col = 'green')
# Plot w/ imprinting
lines(0:99, colSums(sim.H)/sum(sim.H), col = 'magenta', lwd = 4)
legend(15, .07, c('Data - single season', 'Data - all seasons', 'Model - no imprinting', 'Model - imprinting'), pch = c(1, 16, NA, NA), col = c('dodgerblue', 'skyblue', 'white', 'magenta'), bty = 'n', lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3))


paste('Imprinting best')
nll(pars = H1.ests.vv.OAS.H$par, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H1, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE)

paste('No imprinting')
nll(pars = H1.ests.vv.OAS.H$par, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H1, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = FALSE)

paste('OAS')
nll(pars = H1.ests.vv.OAS$par, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H1, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = FALSE)









load('H3ests_vv_OAS_H_2017-09-04.RData')
sim.H = simulate.incidence(alpha.hat = H3.ests.vv.OAS.H$par['alpha.hat'], AA.hat = H3.ests.vv.OAS.H$par['AA.hat'], R0.hat = H3.ests.vv.OAS.H$par['R0.hat'], vv.hat = H3.ests.vv.OAS.H$par['vv.hat'], tau.hat = H3.ests.vv.OAS.H$par['tau.hat'], Hm = H3.ests.vv.OAS.H$par['Hm'], demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE,imprinting = TRUE, wm.in = wm.H3)

load('H3ests_vv_OAS2017-09-03.RData')
sim.OAS = simulate.incidence(alpha.hat = H3.ests.vv.OAS$par['alpha.hat'], AA.hat = H3.ests.vv.OAS$par['AA.hat'], R0.hat = H3.ests.vv.OAS$par['R0.hat'], vv.hat = H3.ests.vv.OAS$par['vv.hat'], tau.hat = H3.ests.vv.OAS$par['tau.hat'], Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H3)

sim.adams.pars = simulate.incidence(alpha.hat = 0.17, AA.hat = 1.20, R0.hat = 2.15, vv.hat = 0.24, tau.hat = 0.93, Hm = NULL, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H3)

# Plot data
par(bg = 'black', fg = 'white', col.main = 'white', col.axis = 'white', col.lab = 'white')
cols = gray.colors(16)[-1]
plot(0:99, H3.inputs[1, ]/sum(H3.inputs[1, ]), ylim = c(0, .08), main = 'H3N2', col = cols[1], ylab = 'Fraction of cases', xlab = 'Case age')
for(ii in 1:nrow(H3.inputs)){ points(0:99, H3.inputs[ii, ]/sum(H3.inputs[ii, ]), col = cols[ii])}
# Plot aggregate data
points(0:99, colSums(H3.inputs)/sum(H3.inputs), pch = 16, cex = .8)
# Plot aggregate OAS
lines(0:99, colSums(sim.OAS)/sum(sim.OAS), col = 'green3', lwd = 2)
lines(0:99, colSums(sim.OAS)/sum(sim.OAS), col = 'greenyellow', lwd = .5)
lines(0:99, colSums(sim.adams.pars)/sum(sim.adams.pars), col = 'blue')
legend(40, .08, 'No imprinting', bty = 'n', lty = 1, col = 'green')
# Plot w/ imprinting
lines(0:99, colSums(sim.H)/sum(sim.H), col = 'dodgerblue3', lwd = 2)
lines(0:99, colSums(sim.H)/sum(sim.H), col = 'lightblue', lwd = .5)
legend(40, .08, c('No imprinting', 'Imprinting'), bty = 'n', lty = 1, col = c('green', 'skyblue'))



# Plot data
cols = 'red'
plot(0:99, H3.inputs[1, ]/sum(H3.inputs[1, ]), col = cols, main = '', ylab = 'Fraction of cases', xlab = 'Case age', ylim = c(0, .09))
for(ii in 1:nrow(H3.inputs)){ points(0:99, H3.inputs[ii, ]/sum(H3.inputs[ii, ]), col = cols)}
mtext('H3N2', line =.5, font = 2)
points(0:99, colSums(H3.inputs)/sum(H3.inputs), pch = 16, col = 'pink2')
# Plot aggregate OAS
lines(0:99, colSums(sim.OAS)/sum(sim.OAS), col = 'white', lwd = 4)
lines(0:99, colSums(sim.OAS)/sum(sim.OAS), col = 'white', lwd = 2.5)
legend(25, .085, c('Data - single season', 'Data - all seasons', 'Model - no imprinting'), pch = c(1, 16, NA), col = c('red', 'pink2', 'white'), bty = 'n', lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3))
#legend(40, .07, 'No imprinting', bty = 'n', lty = 1, col = 'green')
# Plot w/ imprinting
lines(0:99, colSums(sim.H)/sum(sim.H), col = 'purple', lwd = 4)
legend(25, .085, c('Data - single season', 'Data - all seasons', 'Model - no imprinting', 'Model - imprinting'), pch = c(1, 16, NA, NA), col = c('red', 'pink2', 'white', 'purple'), bty = 'n', lty = c(NA, NA, 1, 1), lwd = c(NA, NA, 3, 3))







paste('Imprinting best')
nll(pars = H3.ests.vv.OAS.H$par, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE)

paste('No imprinting')
nll(pars = H3.ests.vv.OAS.H$par, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = FALSE)

paste('Full imprinting')
pars.in = H3.ests.vv.OAS.H$par; pars.in['Hm'] = .5
nll(pars = pars.in, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE)

paste('OAS')
nll(pars = H3.ests.vv.OAS$par, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = FALSE)

