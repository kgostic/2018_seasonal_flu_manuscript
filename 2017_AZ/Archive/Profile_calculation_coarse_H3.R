rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')


#optim(par = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7), nll_prof, Hm = 1, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H3, control = list(maxit = 1000))



Sys.time()
print('Start H3N2 prof valculation OAS Ht')

Hm.vals.H3 = seq(.4, 1, by = .1)
# start a parallel cluster
library(parallel)
no_cores <- min(length(Hm.vals.H3), detectCores()-3)
cl <- makeCluster(no_cores)
# ## export relevant variables
clusterExport(cl, c('Hm.vals.H3', 'epi_final_size', 'simulate.incidence', 'nll_prof', 'demog', 'H3.clusters', 'H3.inputs', 'fraction.vaccinated', 'wm.H3', 'waifw.disaggregate', 'get_first_exposure_probs', 'get_last_exposure_probs'))
#
Hm.vals.H3 = c(.5, 1)
#
# calc profs
H3.OAS.Ht.profs = parLapply(cl = cl, Hm.vals.H3,
                            fun = function(pp.val) optim(par = c(alpha.hat = .08, AA.hat = .89, R0.hat = 3, vv.hat = .5, tau.hat = .7), nll_prof, Hm = pp.val, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H3, control = list(maxit = 1000)) )
stopCluster(cl)

H3.OAS.HT.prof.vals = sapply(H3.OAS.Ht.profs, FUN = function(xx) xx$value)
save(H3.OAS.Ht.profs, H3.OAS.HT.prof.vals, Hm.vals.H3, file = paste('H3_Profs_coarse', Sys.Date(), ".RData", sep = '' ))


