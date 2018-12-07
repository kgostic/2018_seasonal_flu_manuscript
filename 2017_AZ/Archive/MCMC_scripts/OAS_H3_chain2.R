source('MCMC_attempt.R')

set.seed(208)
OAS_H3_chain2 = run.MCMC(data.in = H3.inputs, clust.in = H3.clusters, st.yr = 1968, end.yr = 2014, wm.in = wm.H3, vs = TRUE, OASs = TRUE, is = FALSE, n.steps = 20000, flname = 'OAS_H3_chain2')

save(OAS_H3_chain2, file = 'OAS_H3_chain2.RData')
