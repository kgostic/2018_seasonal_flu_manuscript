## Estimate best values for different models
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')




############################################################
#_________________ Cross-immunity + OAS  ___________________#
#H1N1
Sys.time()
H1.ests.novv.OAS = optim(par = c(alpha.hat = .04, AA.hat = .6, R0.hat = 2.51, vv.hat = NULL, tau.hat = .7), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = FALSE, OAS.switch = TRUE, control = list(maxit = 1500))
beepr::beep(); H1.ests.novv.OAS
save(H1.ests.novv.OAS, file = paste('H1ests_novv_OAS_', Sys.Date(), ".RData", sep = '' ))
#load('H1ests_2017-08-15.RData')
#
#
#
# H3N2
Sys.time()
H3.ests.novv.OAS = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3, vv.hat = NULL, tau.hat = .7), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = FALSE, OAS.switch = TRUE, control = list(maxit = 1500))
beepr::beep(); H3.ests.novv.OAS
save(H3.ests.novv.OAS, file = paste('H3ests_novv_OAS', Sys.Date(), ".RData", sep = '' ))

#load('H3ests_2017-08-02.RData')
############################################################




# 
# 
# ############################################################
# #_________________ Cross-immunity  + OAS + Imprinting  ___________________#
# #H1N1
Sys.time()
H1.ests.novv.OAS.H = optim(par = c(alpha.hat = .04, AA.hat = .6, R0.hat = 2.51, vv.hat = NULL, tau.hat = .7, Hm = .5), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = FALSE, OAS.switch = TRUE, control = list(maxit = 2000), wm.in = wm.H1, imprinting.switch = TRUE)
beepr::beep(); H1.ests.novv.OAS.H
save(H1.ests.novv.OAS.H, file = paste('H1ests_novv_OAS_H_', Sys.Date(), ".RData", sep = '' ))
# #load('H1ests_2017-08-15.RData')
# #
# #
# #
# # H3N2
# Sys.time()
H3.ests.novv.OAS.H = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3, vv.hat = NULL, tau.hat = .7, Hm = .5), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = FALSE, OAS.switch = TRUE, wm.in = wm.H3, imprinting.switch = TRUE, control = list(maxit = 1500))
beepr::beep(); H3.ests.novv.OAS.H
save(H3.ests.novv.OAS.H, file = paste('H3ests_novv_OAS_H_', Sys.Date(), ".RData", sep = '' ))
# 
# #load('H3ests_2017-08-02.RData')
# ############################################################
# 
# 
# 

