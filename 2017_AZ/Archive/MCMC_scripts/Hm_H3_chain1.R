source('MCMC_attempt.R')

Hm_H3_chain1 = run.MCMC(data.in = H3.inputs, clust.in = H3.clusters, st.yr = 1968, end.yr = 2014, wm.in = wm.H3, vs = TRUE, OASs = TRUE, is = TRUE, n.steps = 20000, flname = 'Hm_H3_chain1')

save(Hm_H3_chain1, file = 'Hm_H3_chain1.RData')
