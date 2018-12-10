rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Baseline.R')
set.seed(109)




## Prepare inputs into JAGS model:
######################################################################
# 1. DATA 

# H1.inputs - Matrix giving the number of cases observed in a given age group (columns), 
#             and in a given year (rows). Only years where over 300 cases were observed are included.
# H3.inputs - Same as above for H3N2 data.
# data.in = H1.inputs
# clust.in = H1.clusters
# WAIFW = waifw.disaggregate
# fv.mat = fraction.vaccinated
# vs = TRUE
# OASs = TRUE
# is = TRUE
# st.yr = 1977
# end.yr = 2014
# wm.in = wm.H1












## Define function


run.MCMC = function(data.in, clust.in, WAIFW = waifw.disaggregate, fv.mat = fraction.vaccinated, st.yr, end.yr, wm.in, vs = TRUE, OASs = TRUE, is = TRUE, n.steps = 10000, flname){

start.time = Sys.time()
######################################################################
# 2. INITIAL CONDITIONS
AA = alpha = R0 = vv = tau = Hm = AA.accept = alpha.accept = R0.accept = vv.accept = tau.accept = Hm.accept = vector('numeric', n.steps) # Initialize empty vectors

# Write a function to randomly draw model inits from distributions other than the priors
generate.model.inits = function(){
  c('AA' = rexp(1, 1), 'alpha' = rexp(1, 1.5), 'R0' = runif(1, 1, 4), 'vv' = runif(1, 0, 1), 'tau' = rexp(1, 1.5), 'Hm' = runif(1, 0, 1))
}
pars.inits = generate.model.inits()





#################################################################################
#################################################################################
#################################################################################
# 4. MCMC
# Estimate A, alpha, R0, v, tau, Hm
#################################################################################
#################################################################################
# Set initial conditions
AA[1] = pars.inits['AA']
alpha[1] = pars.inits['alpha']
R0[1] = pars.inits['R0']
vv[1] = pars.inits['vv']
tau[1] = pars.inits['tau']
Hm[1] = pars.inits['Hm']

# nll(pars = c('alpha.hat' = alpha[ss-1], 'AA.hat' = AA[ss-1], 'R0.hat' = R0[ss-1], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, wm.in = wm.H1, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE)



# Set tuning parameters
tune.AA = 1
tune.alpha = 1
tune.R0 = .5
tune.tau = .7
tune.Hm = .1

# Set upper proposal bounds
# Note - 2*tuning must be less than the upper bound
upper.Hm = 1
if(tune.Hm*2 > upper.Hm) warning("2*tuning must be less than the upper Hm bound")

##################         START FOR LOOP HERE      ###########################
for(ss in 2:n.steps){
  
  ## 4.1 MH step for AA
  #      Prior = exponential, lambda = .9
  #      xx = 0:10; plot(xx, dexp(xx, rate = .9))
  AA.new = runif(1, min = -tune.AA, max = tune.AA)+AA[ss-1]#; AA.new # Propose new value
  ## reflect across 0
  AA.new = ifelse(AA.new < 0, -AA.new, AA.new)#; AA.new
  ## Calculate the full conditional at the current value
  log.f.AA.old = -nll(pars = c('alpha.hat' = alpha[ss-1], 'AA.hat' = AA[ss-1], 'R0.hat' = R0[ss-1], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, wm.in = wm.in, vaccinate.switch = vs, OAS.switch = OASs, imprinting.switch = is)+dexp(x = AA[ss-1], rate = .9, log = TRUE)#; log.f.AA.old
  
  
  ## Calculate full conditional at the new proposed value
  log.f.AA.new = -nll(pars = c('alpha.hat' = alpha[ss-1], 'AA.hat' = AA.new,  'R0.hat' = R0[ss-1], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = AA.new, rate = .9, log = TRUE)#; log.f.AA.new
  
  ## Decide whether or not to accept the proposal
  AA.accept[ss] = runif(1) <= min(1, exp(log.f.AA.new-log.f.AA.old))#; AA.accept[ss]
  ## Store the R0 value for this time step
  AA[ss] = ifelse(AA.accept[ss], AA.new, AA[ss-1])#; AA[ss]
  
  
  ## 4.2 MH step for alpha
  #      Prior = exponential, lambda = 1.5
  #      xx = 0:10; plot(xx, dexp(xx, rate = 1.5))
  alpha.new = runif(1, min = -tune.alpha, max = tune.alpha)+alpha[ss-1]#; alpha.new # Propose new value
  ## reflect across 0
  alpha.new = ifelse(alpha.new < 0, -alpha.new, alpha.new)#; alpha.new
  ## Calculate the full conditional at the current value
  log.f.alpha.old = -nll(pars = c('alpha.hat' = alpha[ss-1], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss-1], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = alpha[ss-1], rate = 1.5, log = TRUE)#; log.f.alpha.old
  
  ## Calculate full conditional at the new proposed value
  log.f.alpha.new = -nll(pars = c('alpha.hat' = alpha.new, 'AA.hat' = AA[ss], 'R0.hat' = R0[ss-1], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = alpha.new, rate = 1.5, log = TRUE)#; log.f.alpha.new
  
  ## Decide whether or not to accept the proposal
  alpha.accept[ss] = runif(1) <= min(1, exp(log.f.alpha.new-log.f.alpha.old))#; alpha.accept[ss]
  ## Store the R0 value for this time step
  alpha[ss] = ifelse(alpha.accept[ss], alpha.new, alpha[ss-1])#; alpha[ss]
  
  
  
  ## 4.3 MH step for R0
  #      Prior = exponential, rate = .5
  #      xx = seq(0, 10, by = 1); plot(xx, dexp(xx, .5))
  R0.new = runif(1, min = -tune.R0, max = tune.R0)+R0[ss-1]#; R0.new # Propose new value
  ## reflect across 1
  R0.new = ifelse(R0.new < 1, -R0.new, R0.new)#; R0.new
  ## Calculate the full conditional at the current value
  log.f.R0.old = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss-1], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = R0[ss-1], rate = 0.5, log = TRUE)#; log.f.R0.old
  
  ## Calculate full conditional at the new proposed value
  log.f.R0.new = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0.new, 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = R0.new, rate = 0.5, log = TRUE)#; log.f.R0.new
  
  ## Decide whether or not to accept the proposal
  R0.accept[ss] = runif(1) <= min(1, exp(log.f.R0.new-log.f.R0.old))#; R0.accept[ss]
  ## Store the R0 value for this time step
  R0[ss] = ifelse(R0.accept[ss], R0.new, R0[ss-1])#; R0[ss]
  
  
  
  
  ## 4.4 MH step for vv
  #      Prior = beta, shape 1 = 1.75, shape 2 = 1.5
  #      xx = seq(0, 1, by = .01); plot(xx, dbeta(xx, shape1 = 1.75, 1.5))
  vv.new = runif(1)#; vv.new # Propose new value
  ## Calculate the full conditional at the current value
  log.f.vv.old = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss], 'vv.hat' = vv[ss-1], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dbeta(vv[ss-1], shape1 = 1.75, shape2 = 1.5, log = TRUE)#; log.f.vv.old
  
  ## Calculate full conditional at the new proposed value
  log.f.vv.new = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss], 'vv.hat' = vv.new, 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dbeta(vv.new, shape1 = 1.75, shape2 = 1.5, log = TRUE)#; log.f.vv.new
  
  ## Decide whether or not to accept the proposal
  vv.accept[ss] = runif(1) <= min(1, exp(log.f.vv.new-log.f.vv.old))#; vv.accept[ss]
  ## Store the vv value for this time step
  vv[ss] = ifelse(vv.accept[ss], vv.new, vv[ss-1])#; vv[ss]
  
  
  
  
  ## 4.5 MH step for tau
  #      Prior = exponential, rate = 1
  #      xx = seq(0, 10, by = 1); plot(xx, dexp(xx, 1))
  tau.new = runif(1, min = -tune.tau, max = tune.tau)+tau[ss-1]#; tau.new # Propose new value
  ## reflect across 0
  tau.new = ifelse(tau.new < 0, -tau.new, tau.new)#; tau.new
  ## Calculate the full conditional at the current value
  log.f.tau.old = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss], 'vv.hat' = vv[ss], 'tau.hat' = tau[ss-1], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = tau[ss-1], rate = 1, log = TRUE)#; log.f.tau.old
  
  ## Calculate full conditional at the new proposed value
  log.f.tau.new = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss], 'vv.hat' = vv[ss], 'tau.hat' = tau.new, 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)+dexp(x = tau.new, rate = 1, log = TRUE)#; log.f.tau.new
  
  ## Decide whether or not to accept the proposal
  tau.accept[ss] = runif(1) <= min(1, exp(log.f.tau.new-log.f.tau.old))#; tau.accept[ss]
  ## Store the tau value for this time step
  tau[ss] = ifelse(tau.accept[ss], tau.new, tau[ss-1])#; tau[ss]
  
 
  if(is == TRUE){ 
  ## 4.6 MH step for Hm
  #      Prior = unif[0, 1],  = 1 for all values
  Hm.new = runif(1, min = -tune.Hm, max = tune.Hm)+Hm[ss-1]#; Hm.new # Propose new value
  ## reflect across 0
  Hm.new = ifelse(Hm.new < 0, -Hm.new, Hm.new)#; Hm.new
  ## reflect across 1
  Hm.new = ifelse(Hm.new > 1, 1-Hm.new, Hm.new)#; Hm.new
  ## Calculate the full conditional at the current value
  log.f.Hm.old = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss], 'vv.hat' = vv[ss], 'tau.hat' = tau[ss], 'Hm' = Hm[ss-1]), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)#; log.f.Hm.old
  
  ## Calculate full conditional at the new proposed value
  log.f.Hm.new = -nll(pars = c('alpha.hat' = alpha[ss], 'AA.hat' = AA[ss], 'R0.hat' = R0[ss], 'vv.hat' = vv[ss], 'tau.hat' = tau[ss], 'Hm' = Hm.new), demog = demog, clusters = clust.in, dat.in = data.in, wm.in = wm.in, st.yr = st.yr, end.yr = end.yr, plot = FALSE, WAIFW.in = WAIFW, vax.matrix = fv.mat, vaccinate = vs, OAS = OASs, imprinting = is)#; log.f.Hm.new
  
  ## Decide whether or not to accept the proposal
  Hm.accept[ss] = runif(1) <= min(1, exp(log.f.Hm.new-log.f.Hm.old))#; Hm.accept[ss]
  ## Store the Hm value for this time step
  Hm[ss] = ifelse(Hm.accept[ss], Hm.new, Hm[ss-1])#; Hm[ss]
  
  }
  
  
  if(ss %% 100 == 0){ # Every 500 steps, save the outputs
  ests = rbind(alpha, AA, R0, vv, tau, Hm); rownames(ests) = c('alpha', 'AA', 'R0', 'vv', 'tau', 'Hm')
  acceptance = rbind(alpha.accept, AA.accept, R0.accept, vv.accept, tau.accept, Hm.accept); rownames(ests) = c('alpha', 'AA', 'R0', 'vv', 'tau', 'Hm')
  save(ests, acceptance, file = paste(flname, '_', ss, 'steps_', Sys.Date(), '.RData', sep = ''))
  print(start.time - Sys.time())
  }
}
########################  End Run  ##########################

# Return
ests = rbind(alpha, AA, R0, vv, tau, Hm); rownames(ests) = c('alpha', 'AA', 'R0', 'vv', 'tau', 'Hm')
acceptance = rbind(alpha.accept, AA.accept, R0.accept, vv.accept, tau.accept, Hm.accept); rownames(ests) = c('alpha', 'AA', 'R0', 'vv', 'tau', 'Hm')

print(start.time - Sys.time())
list(ests = ests, accept = acceptance)
}


# 
# 
#run.MCMC(data.in = H1.inputs, clust.in = H1.clusters, st.yr = 1977, end.yr = 2014, wm.in = wm.H1, vs = TRUE, OASs = TRUE, is = TRUE, n.steps = 10, flname = 'Hm_H1_chain1')
# 
# run.MCMC(data.in = H1.inputs, clust.in = H1.clusters, st.yr = 1977, end.yr = 2014, wm.in = wm.H1, vs = TRUE, OASs = TRUE, is = FALSE, n.steps = 10, flname = 'OAS_H1_chain1')
# 
# run.MCMC(data.in = H3.inputs, clust.in = H3.clusters, st.yr = 1977, end.yr = 2014, wm.in = wm.H3, vs = TRUE, OASs = TRUE, is = TRUE, n.steps = 10, flname = 'Hm_H3_chain1')
# 
# run.MCMC(data.in = H3.inputs, clust.in = H3.clusters, st.yr = 1977, end.yr = 2014, wm.in = wm.H3, vs = TRUE, OASs = TRUE, is = FALSE, n.steps = 10, flname = 'OAS_H3_chain1')
# 
