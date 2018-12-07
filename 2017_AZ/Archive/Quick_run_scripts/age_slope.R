## Estimate best values for different models
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')

### age ONLY
load('H3ests_vv_OAS2017-09-03.RData')






### age and imprinting combined, with sloped age pattern
Sys.time()
bb = H3.ests.vv.OAS$par; bb['alpha.hat'] = 1.5
H3.ests.vv.OAS.age_slope = optim(par = c(H3.ests.vv.OAS$par, age.par = 20), fn = nll_age_imp, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), slope = TRUE, age.switch = TRUE, imp.switch = FALSE)
beepr::beep(); H3.ests.vv.OAS.age_slope
save(H3.ests.vv.OAS.age_slope, file = paste('H3ests_vv_OAS_age_slope', Sys.Date(), ".RData", sep = '' ))

