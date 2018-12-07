source('Profile_functions.R')

print('Start H1N1 prof OAS Imprinting Transmission')
Sys.time()
#H1N1

# start a parallel cluster
library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
# ## export relevant variables
clusterExport(cl, c('epi_final_size', 'simulate.incidence', 'nll_prof', 'demog', 'H1.clusters', 'H1.inputs', 'fraction.vaccinated', 'wm.H1', 'waifw.disaggregate', 'get_first_exposure_probs', 'get_last_exposure_probs'))
#
Hm.vals.H1 = seq(0, 1, length = 3*no_cores)
#
# calc profs
H1.OAS.Ht.profs = parLapply(cl = cl, Hm.vals.H1,
                            fun = function(pp.val) optim(par = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7), nll_prof, Hm = pp.val, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H1, control = list(maxit = 1000)) )

stopCluster(cl)

H1.OAS.HT.prof.vals = sapply(H1.OAS.Ht.profs, FUN = function(xx) xx$value)
save(H1.OAS.Ht.profs, H1.OAS.HT.prof.vals, Hm.vals.H1, file = paste('H1_Profs_', Sys.Date(), ".RData", sep = '' ))





Sys.time()
print('Start H3N2 prof valculation OAS Ht')
# start a parallel cluster
library(parallel)
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
# ## export relevant variables
clusterExport(cl, c('epi_final_size', 'simulate.incidence', 'nll_prof', 'demog', 'H3.clusters', 'H3.inputs', 'fraction.vaccinated', 'wm.H3', 'waifw.disaggregate', 'get_first_exposure_probs', 'get_last_exposure_probs'))
#
Hm.vals.H3 = seq(0, 1, length = 3*no_cores)
#
# calc profs
H3.OAS.Ht.profs = parLapply(cl = cl, Hm.vals.H3,
                            fun = function(pp.val) optim(par = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7), nll_prof, Hm = pp.val, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H3, control = list(maxit = 1000)) )
stopCluster(cl)

H3.OAS.HT.prof.vals = sapply(H3.OAS.Ht.profs, FUN = function(xx) xx$value)
save(H3.OAS.Ht.profs, H3.OAS.HT.prof.vals, Hm.vals.H3, file = paste('H3_Profs_', Sys.Date(), ".RData", sep = '' ))










