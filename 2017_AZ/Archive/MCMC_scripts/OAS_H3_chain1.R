source('MCMC_attempt.R')

OAS_H3_chain1 = run.MCMC(data.in = H3.inputs, clust.in = H3.clusters, st.yr = 1968, end.yr = 2014, wm.in = wm.H3, vs = TRUE, OASs = TRUE, is = FALSE, n.steps = 20000, flname = 'OAS_H3_chain1')

save(OAS_H3_chain1, file = 'OAS_H3_chain1.RData')
