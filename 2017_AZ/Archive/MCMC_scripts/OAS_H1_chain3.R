source('MCMC_attempt.R')

set.seed(32)
OAS_H1_chain3 = run.MCMC(data.in = H1.inputs, clust.in = H1.clusters, st.yr = 1977, end.yr = 2014, wm.in = wm.H1, vs = TRUE, OASs = TRUE, is = FALSE, n.steps = 20000, flname = 'OAS_H1_chain3')

save(OAS_H1_chain3, file = 'OAS_H1_chain3.RData')
