## Estimate best values for different models
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')

# ############################################################
# #_________________ Cross-immunity only ___________________#
# #H1N1
# Sys.time()
# H1.ests = optim(par = c(alpha.hat = .04, AA.hat = .6, R0.hat = 2.51), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = FALSE, OAS.switch = FALSE, imprinting.switch = FALSE)
# beepr::beep(); H1.ests
# save(H1.ests, file = paste('H1ests_', Sys.Date(), ".RData", sep = '' ))
# #load('H1ests_2017-08-15.RData')
# #
# #
# #
# # H3N2
# Sys.time()
# H3.ests = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = FALSE, OAS.switch = FALSE)
# beepr::beep(); H3.ests
# save(H3.ests, file = paste('H3ests_', Sys.Date(), ".RData", sep = '' ))
# 
# #load('H3ests_2017-08-02.RData')
# ############################################################
# 
# 
# 
# 
# 
# 
# 
# ############################################################
# #_________________ Cross-immunity + vaccination  ___________________#
# #H1N1
# Sys.time()
# H1.ests.vv = optim(par = c(alpha.hat = .04, AA.hat = .6, R0.hat = 2.51, vv.hat = .5), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = FALSE, control = list(maxit = 1500))
# beepr::beep(); H1.ests.vv
# save(H1.ests.vv, file = paste('H1ests_vv_', Sys.Date(), ".RData", sep = '' ))
# #load('H1ests_2017-08-15.RData')
# #
# #
# #
# # H3N2
# Sys.time()
# H3.ests.vv = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3, vv.hat = .5), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = FALSE, control = list(maxit = 1500))
# beepr::beep(); H3.ests.vv
# save(H3.ests.vv, file = paste('H3ests_vv_', Sys.Date(), ".RData", sep = '' ))
# 
# #load('H3ests_2017-08-02.RData')
# ############################################################


############################################################
#_________________ Cross-immunity + vaccination + OAS  ___________________#
#H1N1
Sys.time()
H1.ests.vv.OAS = optim(par = c(alpha.hat = .04, AA.hat = .6, R0.hat = 2.51, vv.hat = .5, tau.hat = .7), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 1500))
beepr::beep(); H1.ests.vv.OAS
save(H1.ests.vv.OAS, file = paste('H1ests_vv_OAS_', Sys.Date(), ".RData", sep = '' ))
#load('H1ests_2017-08-15.RData')
#
#
#
# H3N2
Sys.time()
H3.ests.vv.OAS = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3, vv.hat = .5, tau.hat = .7), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2015, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 1500))
beepr::beep(); H3.ests.vv.OAS
save(H3.ests.vv.OAS, file = paste('H3ests_vv_OAS', Sys.Date(), ".RData", sep = '' ))

#load('H3ests_2017-08-02.RData')
############################################################




# 
# 
# ############################################################
# #_________________ Cross-immunity + vaccination + OAS + Imprinting  ___________________#
# #H1N1
Sys.time()
load('H1ests_vv_OAS_2017-09-03.RData')
pp = H1.ests.vv.OAS$par
H1.ests.vv.OAS.H = optim(par = c(pp, Hm = 1), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), wm.in = wm.H1, imprinting.switch = TRUE)
beepr::beep(); H1.ests.vv.OAS.H
save(H1.ests.vv.OAS.H, file = paste('H1ests_vv_OAS_H_', Sys.Date(), ".RData", sep = '' ))
# #load('H1ests_2017-08-15.RData')
# #
# #
# #
# # H3N2
# Sys.time()
H3.ests.vv.OAS.H = optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3, vv.hat = .5, tau.hat = .7, Hm = .5), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2015, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, wm.in = wm.H3, imprinting.switch = TRUE, control = list(maxit = 1500))
beepr::beep(); H3.ests.vv.OAS.H
save(H3.ests.vv.OAS.H, file = paste('H3ests_vv_OAS_H_', Sys.Date(), ".RData", sep = '' ))
# 
# #load('H3ests_2017-08-02.RData')
# ############################################################
# 
# 










##### CASE ASCERTAINMENT MODELS!!!

## See quick plots folder to run these!

### age and imprinting combined
load('H3ests_vv_OAS2017-09-03.RData')
Sys.time()
H3.ests.vv.OAS.age.imp = optim(par = c(H3.ests.vv.OAS$par, HH = 3, age.par = 10), fn = nll_age_imp, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000))
beepr::beep(); H3.ests.vv.OAS.age.imp
save(H3.ests.vv.OAS.age.imp, file = paste('H3ests_vv_OAS_age_imp', Sys.Date(), ".RData", sep = '' ))


### age ONLY
Sys.time()
H3.ests.vv.OAS.age = optim(par = c(H3.ests.vv.OAS$par, age.par = 10), fn = nll_age_imp, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2015, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), age.switch = TRUE, imp.switch = FALSE)
beepr::beep(); H3.ests.vv.OAS.H.age
save(H3.ests.vv.OAS.H.age, file = paste('H3ests_vv_OAS_H_age', Sys.Date(), ".RData", sep = '' ))




### imp only
Sys.time()
H3.ests.vv.OAS.imp = optim(par = c(H3.ests.vv.OAS$par, HH = 3), fn = nll_age_imp, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), age.switch = FALSE, imp.switch = TRUE)
beepr::beep(); H3.ests.vv.OAS.imp
save(H3.ests.vv.OAS.imp, file = paste('H3ests_vv_OAS_imp', Sys.Date(), ".RData", sep = '' ))




### age and imprinting combined, with sloped age pattern
Sys.time()
H3.ests.vv.OAS.age_slope.imp = optim(par = c(H3.ests.vv.OAS$par, HH = 3, age.par = 10), fn = nll_age_imp, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), slope = TRUE, age.switch = TRUE, imp.switch = TRUE)
beepr::beep(); H3.ests.vv.OAS.H.age_slope.imp
save(H3.ests.vv.OAS.H.age_slope.imp, file = paste('H3ests_vv_OAS_H_age_slope_imp', Sys.Date(), ".RData", sep = '' ))

