rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')

print('Start H1N1 prof OAS Imprinting Transmission')
Sys.time()
#H1N1

# start a parallel cluster
Hm.vals.H1 = seq(.4, 1, by = .1)
library(parallel)
no_cores <-  min(length(Hm.vals.H1), detectCores()-3)
cl <- makeCluster(no_cores)
# ## export relevant variables
clusterExport(cl, c('Hm.vals.H1', 'epi_final_size', 'simulate.incidence', 'nll_prof', 'demog', 'H1.clusters', 'H1.inputs', 'fraction.vaccinated', 'wm.H1', 'waifw.disaggregate', 'get_first_exposure_probs', 'get_last_exposure_probs'))
#

#
# calc profs
H1.OAS.Ht.profs = parLapply(cl = cl, Hm.vals.H1,
                            fun = function(pp.val) optim(par =c(alpha.hat = .04, AA.hat = .6, R0.hat = 2.51, vv.hat = .5, tau.hat = .7), nll_prof, Hm = pp.val, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H1, control = list(maxit = 1000)) )



stopCluster(cl)

H1.OAS.HT.prof.vals = sapply(H1.OAS.Ht.profs, FUN = function(xx) xx$value)
save(H1.OAS.Ht.profs, H1.OAS.HT.prof.vals, Hm.vals.H1, file = paste('H1_Profs_coarse', Sys.Date(), ".RData", sep = '' ))



