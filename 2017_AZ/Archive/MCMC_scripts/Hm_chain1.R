source('MCMC_attempt.R')

Hm_H1_chain1 = run.MCMC(data.in = H1.inputs, clust.in = H1.clusters, st.yr = 1977, end.yr = 2014, wm.in = wm.H1, vs = TRUE, OASs = TRUE, is = TRUE, n.steps = 20000, flname = 'Hm_H1_chain1')

save(Hm_H1_chain1, file = 'Hm_H1_chain1.RData')
