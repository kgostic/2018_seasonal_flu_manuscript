
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Profile_functions.R')

H3.5 = optim(par = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7), nll_prof, Hm = .5, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H3, control = list(maxit = 1000))

save(H3.5, file = 'H3_prof_5.RData')