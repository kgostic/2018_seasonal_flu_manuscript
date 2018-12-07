source('MCMC_attempt.R')

set.seed(32)
Hm_H3_chain3 = run.MCMC(data.in = H3.inputs, clust.in = H3.clusters, st.yr = 1968, end.yr = 2014, wm.in = wm.H3, vs = TRUE, OASs = TRUE, is = TRUE, n.steps = 20000, flname = 'Hm_H3_chain3')

save(Hm_H3_chain3, file = 'Hm_H3_chain3.RData')
