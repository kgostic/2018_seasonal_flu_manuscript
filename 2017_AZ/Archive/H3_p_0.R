
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')
rm(wm.H1)

H3_0 = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 10, vv.hat = .5, tau.hat = .7), nll_prof, Hm = 0, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H3, control = list(maxit = 1000), penalty.in = 2500)

save(H3_0, file = 'H3_prof_0.RData')