## Estimate best values for different models
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')

### age ONLY
load('H3ests_vv_OAS2017-09-03.RData')







### imp only
Sys.time()
H3.ests.vv.OAS.imp = optim(par = c(H3.ests.vv.OAS$par, HH = 3), fn = nll_age_imp, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = wm.H3, vaccinate.switch = TRUE, OAS.switch = TRUE, control = list(maxit = 2000), age.switch = FALSE, imp.switch = TRUE)
beepr::beep(); H3.ests.vv.OAS.imp
save(H3.ests.vv.OAS.imp, file = paste('H3ests_vv_OAS_imp', Sys.Date(), ".RData", sep = '' ))

