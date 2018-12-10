rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')

############################################################
#_________________ Cross-immunity + vaccination + OAS + Imprinting  ___________________#
#H1N1
# Sys.time()
# H1.ests.vv.OAS.H = optim(par = c(alpha.hat = .01, AA.hat = .74, R0.hat = 3.24, vv.hat = .01, tau.hat = .3, Hm = .9), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), wm.in = wm.H1, imprinting.switch = TRUE)
# beepr::beep()
# H1.ests.vv.OAS.H
# save(H1.ests.vv.OAS.H, file = paste('H1ests_vv_OAS_H_', Sys.Date(), ".RData", sep = '' ))
# #load('H1ests_2017-08-15.RData')
#
#
#
# H3N2
Sys.time()
H3.ests.vv.OAS.H = optim(par = c(alpha.hat = .09, AA.hat = .72, R0.hat = 1.1, vv.hat = .0001, tau.hat = 1.7, Hm = .9), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, wm.in = wm.H3, imprinting.switch = TRUE, control = list(maxit = 1500))
beepr::beep()
H3.ests.vv.OAS.H
save(H3.ests.vv.OAS.H, file = paste('H3ests_vv_OAS_H_', Sys.Date(), ".RData", sep = '' ))

#load('H3ests_2017-08-02.RData')
############################################################



