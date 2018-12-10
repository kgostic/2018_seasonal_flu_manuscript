rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Kucharski_SIRfinalsize/R/epi_final_size.R')

#### Simulate the age-specific incidence of seasonal influenza
#### Begin simulation in 1968 for H3N2 and in 1977 for H1N1
#### Track the proportion of each age group with past exposure to all clusters {1, 2, ... n} in the set Y
#### Assume all seasonal variants in the same cluster are antigenically identical, and decreasing cross immunity between strains in increasingly distant clusters.
#### Use model formulation from http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002741
#### Assume cluster jumps as in http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s005&type=supplementary



#_________ SET INPUT PARAMETERS for H1N1  __________# 
# Source: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s006&type=supplementary
# R0.hat = 1.2
# AA.hat = 0.34
# alpha.hat = 0.25
# nu.hat = 0.06
#####################################################




#_________ Define demographic vectors __________# 
# First, for each individual age from 0-99
#demog = c(rep(1/80.5, 60), -1/3220*((60:99)-60)+1/80.5); sum(demog) # This is a rough demographic curve 

# Import demographic data from the USA for all available years (1980-2020)
demog.raw = read.csv('USA_demography_1980_2020_formatted.csv', header = T)
# Drop individuals over age 99
demog = demog.raw[,-101]/rowSums(demog.raw[,-101]) # Normalize
rownames(demog) = 1980:2020; colnames(demog) = 0:99


# Then, for each 5-year age group in POLYMOD
age.grp.index = cbind('lower' = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)+1, 'upper' = c(4, 9, 14, 19, 24, 29, 34, 39, 44, 49, 54, 59, 64, 69, 99)+1)
d.group.collapsed = matrix(NA, nrow = nrow(demog), ncol = nrow(age.grp.index)) 
for(ii in 1:nrow(demog)){
  d.group.collapsed[ii, ] = apply(age.grp.index, MARGIN = 1, function(x) sum(demog[ii, x[1]:x[2]])) # Fraction of population in each age group
}
# Repeat each sum once for each age in the group
d.group = t(apply(d.group.collapsed, 1, function(x) rep(x, c(rep(5, 14), 30))))
# Below, use d.group to rescale the contact matrix for each year 



#_________ Define contact (WAIFW) matrix __________# 
# Source: https://doi.org/10.1371/journal.pmed.0050074.st005
waifw.aggregate = as.matrix(read.csv('WAIFW_GB_Physical.csv', header = TRUE))
# Each column is the age group initiating contact, each row is the age group of the contact made
# Expand so that rates correspond to single years and not five-year bins
waifw.raw = (matrix(rep(waifw.aggregate,  each = 5), nrow = nrow(waifw.aggregate), ncol = ncol(waifw.aggregate)*5, byrow = T)) # Expand acroww columns
waifw.raw = t(matrix(rep(waifw.raw,  each = 5), nrow = nrow(waifw.raw)*5, ncol = ncol(waifw.raw))) # Expand across rows
colnames(waifw.raw) = paste('from', 0:99, sep = '.')
rownames(waifw.raw) = paste('to', 0:99, sep = '.')


### Disaggregate the transmission matrix from 5-year age groups to single-year age groups
# Approach - let a and b represent 5-year age groups and i, j represent single ages falling within a and b, respectively. Then:
# m_ab*Pj/Pb = m_aj   Multiply by the faction of contacts with group b that belong to age j
# m_aj*Pi/Pa = m_ij   Multiply by the faction of contacts from group a that initiated by individuals in age i
waifw.disaggregate = array(NA, dim = c(nrow(waifw.raw), ncol(waifw.raw), nrow(demog)), dimnames = list(NULL, NULL, 1980:2020))
# Scale m_ij using the demographic vectors specific to each year
for(yy in 1:nrow(demog)){
  waifw.disaggregate[,,yy] = t( t(waifw.raw*as.numeric(demog[yy,]/d.group[yy,]))*as.numeric(demog[yy,]/d.group[yy,]) )
}
image(0:99, 0:99, log(as.matrix(waifw.disaggregate[,,10])), main = 'Disaggregate waifw')
# Dimensions of m_ij:
#  * Rows represent age of contact source
#  * Columns represent age of contact target
#  * Dimension 3 represents year (because demographic disaggregation is different based on each year's demographic curve)



#_________ Define clusters __________#
# Source 1: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s005&type=supplementary
# Source 2: WHO reports, see Cluster_transition_notes.xlsx
H1.clusters = data.frame('year' = 1977:2017, 
                         'jump.year' = c(1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,1,0,0),
                         'cluster' = c(1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11))

H3.clusters = data.frame('year' = 1968:2017, 
                         'jump.year' = c(1,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,0,0,1,1,1,0,1,0,1,0,1,0,1,0,0,0,0),
                         'cluster' = c(1,1,1,1,2,2,2,2,3,3,4,4,5,5,5,5,5,5,5,6,6,7,7,7,8,8,8,8,9,10,10,10,10,10,10,11,12,13,13,14,14,15,15,16,16,17,17,17,17,17))



#_________ Define a funciton to calculate the prob of transmission   __________#
#_________       transmission given past exposure to a given set     __________#

# sets is a binary vector indicating which clusters are contained in the set
# challenge.value gives the index of the challenge strain
get.p.tns = function(set, challenge.value, AA.hat, alpha.hat){
  if(any(set == 1)){ # If any prior exposures
    distances = abs(challenge.value - which(set == 1))
    max(0, min(1-AA.hat*exp(-alpha.hat*distances)) )
  }else{ # Else, no prior exposures and fully transmisssible
    1
  }
}

#____________________ Load data on USA vaccination ______________________#
raw.frac.vaccinated = read.csv('USA_vaccination_coverage.csv', skip = 1, header = TRUE)
fraction.vaccinated = matrix( c(raw.frac.vaccinated[,2], # age 0
                                rep(raw.frac.vaccinated[,3], 17), # age 1-17
                                rep(raw.frac.vaccinated[,4], 32), # age 18-49 
                                rep(raw.frac.vaccinated[,5], 15), # age 50-64
                                rep(raw.frac.vaccinated[,6], 35)),# age 65-99
                              nrow = length(1968:2017), ncol = 100, dimnames = list(1968:2017, 0:99))/100








############################################################
#____________________ Load NCBI data ______________________#
#
############################################################
H1_NCBI = read.csv('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/NCBI_H1_Data.csv', stringsAsFactors = FALSE)
H3_NCBI = read.csv('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/NCBI_H3_Data.csv', stringsAsFactors = FALSE)

## Load Arizona data
H1_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H1.csv')
H3_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H3.csv')


################# 
reformat.data = function(dat.in){
  yrs = unique(dat.in$Year)
  ages = 0:99
  out.mat = matrix(NA, nrow = length(yrs), ncol = length(ages), dimnames = list(yrs, ages))
  for(yy in 1:length(yrs)){
    for(aa in 1:length(ages)){
      out.mat[yy, aa] = sum(dat.in$Year == yrs[yy] & dat.in$Age == ages[aa])
    }
  }
  out.mat
}

## Reformatted data inputs
H1.inputs = reformat.data(H1_AZ)
H3.inputs = reformat.data(H3_AZ)









#____________________ Get imprinting weights ______________________#
# source('Infection.age.structure.vaccination.R')
# vaccination.matrix.raw = read.csv('~/Dropbox/R/ImmuneAgeStructure/ForSimulations/Demography_files/vaccinationRates_ages0to9.csv')
# vaccination.matrix = cbind(year = vaccination.matrix.raw$year, vaccination.matrix.raw[,-1]*0.6) #Assume 0.6 in naiive children
# wts = get.type.weights.AB.vaccination_2(years.out = 1968:2017, Countries.out = c('USA'), type = 5, vax.matrix = vaccination.matrix, earliest.birth.year = 1869)
# 
# w1 = wts$weights.master.1 
# w2 = wts$weights.master.2
# w3 = wts$weights.master.3
# wn = wts$weights.master.naiive
# 
# ## Reformat imprinting data
# wm.H1 = wo.H1 = fraction.vaccinated*0 # Initalize
# for(ii in 1:nrow(fraction.vaccinated)){
#   yr = as.numeric(rownames(fraction.vaccinated))[ii]
#   wm.H1[ii, ] = as.numeric((w1+w2)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
#   wo.H1[ii, ] = as.numeric((w3+wn)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
# }
# wm.H3 = wo.H3 = fraction.vaccinated*0 # Initalize
# for(ii in 1:nrow(fraction.vaccinated)){
#   yr = as.numeric(rownames(fraction.vaccinated))[ii]
#   wm.H3[ii, ] = as.numeric((w3)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
#   wo.H3[ii, ] = as.numeric((w1+w2+wn)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
# }
# 
# save(wm.H1, wm.H3, wo.H1, wo.H3, file = 'Seasonal_group_weights.RData')
load('Seasonal_group_weights.RData')
############################################################
#######              END LOADING DATA                 ######
############################################################


# For function debugging
# alpha.hat = .4
# AA.hat = 3
# R0.hat = 3
# tau.hat = .4
# vv.hat = .5
# Hm = .5
# wm.in = as.matrix(wm.H1)
# wo.in = as.matrix(wo.H1)
# start.year = 1977
# end.year = 2014
# clusters = H1.clusters




############################################################
#____________ Simulation function  ______________#
#
############################################################
simulate.incidence = function(alpha.hat, AA.hat, R0.hat, vv.hat, tau.hat, Hm, demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, fraction.vaccinated, wm.in, wo.in){
  
  
  # Use epi_final_size function
  # Initialize at 1977 and itterate forward to 2014 (the last year where we have data)
  years = start.year:end.year
  wm.use = wm.in[as.character(years), ]
  wo.use = wo.in[as.character(years), ]
  
  p.infected = matrix(NA, nrow = length(years), ncol = 100, dimnames = list(years, 0:99)) # Declare output storage
  fraction.vaccinated = fraction.vaccinated[as.character(years), ]
  
  ## Initialize a matrix to keep track of what fraction of each age group has previous exposure to each strain
  frac.prior.exposure = matrix(0, nrow = max(clusters$cluster), ncol = 100); colnames(frac.prior.exposure) = paste('age.', 0:99, sep = ''); rownames(frac.prior.exposure) = paste('cluster.', 1:max(clusters$cluster), sep = '')
  ## Rows define ages, columns define sets and entries the define the fraction of the population in age group a with previous exposure to strain X
  
  
  # These matrices track which strains belong in each set
  # Each row is a different possible set
  # Columns take values 0 or 1 to indicate whether the set contains a strain from each cluster
  YYs = matrix(0, nrow = 1, ncol = max(clusters$cluster), dimnames = list(NULL, paste('cluster', 1:max(clusters$cluster))))
  # First row defines the null set
  # Subsequent rows will be added within the simulation
  
  
  #_________ Caculate final sizes __________#
  # 1977 Step
  # Reduce the initial fraction based on childhood imprinting
  frac.susceptible.before.imprinting.and.vaccination = rep(1, 100) # Initially, assume everyone is susecptible
  frac.susceptible.before.vaccination = frac.susceptible.before.imprinting.and.vaccination * (Hm*wm.use[1, ]+wo.use[1, ])
  
  # And also based on vaccination
  frac.susceptible = frac.susceptible.before.vaccination * (1 - fraction.vaccinated[1,] * vv.hat)
  p.infected[1, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = frac.susceptible)
  # Optional plot of results
  #plot(0:99, p.infected[1,]*demog[1,], main = years[1], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
  
  # Update the sets of possible past exposures
  YYs = rbind(YYs, YYs)
  # Then add exposure to the new strain to each of the new sets
  YYs[2, 1] = 1
  # Update the fraction of people with past exposure to the strain of interest
  frac.prior.exposure[1, ] = p.infected[1,]
  # Account for aging
  frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
  frac.prior.exposure[,1] = 0
  
  ## Assume no OAS effects in the first year, as no one has imprinted yet
  
  
  
  # 1978 - 2017 epidemics
  for(yy in 2:length(years)){
    past.clusters = unique(clusters$cluster[1:yy-1]) # Vector of all clusters that have circulated previously
    n.past.clusters = length(past.clusters)
    n.past.sets = 2^n.past.clusters
    challenge.strain = clusters$cluster[yy]
    
    #___________ Calculate the fraction susceptible _________#
    #S_Ys is a matrix that gives the fraction of each age group with past exposure to set Y
    S_Ys = matrix(NA, n.past.sets, 100)
    
    for(ss in 1:nrow(YYs)){
      ## Multiply across all probabilities of having been previously exposed to each strain in the set, and of having not been previously exposed to each strain not in the set
      #S_Ys[ss, ] = apply((YYs[ss, ]*frac.prior.exposure+(1-YYs[ss, ])*(1-frac.prior.exposure)), 2, prod)
      S_Ys[ss, ] = exp(colSums(log( YYs[ss, ]*frac.prior.exposure+(1-YYs[ss, ])*(1-frac.prior.exposure) )))
    }
    
    #if(any(round(S_Ys, 6) != round(S_Ys.2, 6))) warning('New calculation not working')
    #if(!any(round(S_Ys, 6) != round(S_Ys.2, 6))) print('New calculation works!')
    
    if(any(round(colSums(S_Ys, na.rm = TRUE), 10) != 1)) warning('Faulty S_Y calculation, fractions do not sum to 1')
    
    
    
    # Get the fraction of each age group that can transmit the challenge strain
    p.tns = apply(YYs, 1, function(ss) get.p.tns(ss, challenge.strain, AA = AA.hat, alpha = alpha.hat))
    
    # Reduce the initial fraction based on childhood imprinting
    frac.susceptible.before.imprinting.and.vaccination = as.vector(t(S_Ys) %*% p.tns) # Gives a vector of the probability that each age group will transmit, given exposure histories
    frac.susceptible.before.vaccination = frac.susceptible.before.imprinting.and.vaccination * (Hm*wm.use[yy, ]+wo.use[yy, ]) # Reduce the susceptible fraction based on prob of imprinting protection
    # Add the effects of vaccination
    frac.susceptible = frac.susceptible.before.vaccination * (1 - fraction.vaccinated[yy,] * vv.hat)
    # vaccinated.fraction is an input vector describing the probability of vaccination in each age group
    # vv is a parameter to be fit, which describes vaccine effectiveness
    
    
    if(years[yy] <= 1980){ # 1980 is the first year in which demographic information was available, so use demog[1, ] vector for all years from 1977:1980
      p.infected[yy, ] = epi_final_size(r0 = R0.hat, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = frac.susceptible)
      # plot(0:99, p.infected[yy,]*demog[yy,], main = years[yy], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
      
    }else{
      rr = 2
      p.infected[yy, ] = epi_final_size(r0 = R0.hat, contact_matrix = waifw.disaggregate[,,rr], demography_vector = demog[rr,], prop_suscep = frac.susceptible)
      #plot(0:99, p.infected[yy,]*demog[yy,], main = years[yy], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
      rr = rr+1
    }
    
    # Model original antigenic sin
    first.cluster = apply(YYs, 1, function(yy) min(which(yy == 1))) # In each possible set of first exposures, which was the first strain to circulate?
    first.cluster[1] = NA # Correct the null set
    # Generate a vector listing the strains which will prevent formation of new immune memory to the challenge strain based on OAS
    OAS.strains = which( apply(diag(nrow = n.past.clusters, ncol = n.past.clusters), 1, function(xx) get.p.tns(set = xx, challenge.value = challenge.strain, AA.hat = AA.hat, alpha.hat = alpha.hat)) < tau.hat )
    # Exclude sets subject to OAS from the overall list of possible exposures
    S_Ys.no.OAS = S_Ys[-which(first.cluster %in% OAS.strains), ]
    # What fraction of each age group is NOT prevented from forming new immunity by OAS:
    if(is.null(nrow(S_Ys.no.OAS))){
    OAS.free.frac = S_Ys.no.OAS
    }else{
    OAS.free.frac = colSums(S_Ys.no.OAS)
    }
    
    # Update the sets of possible past exposures
    if(challenge.strain > max( which(colSums(YYs) > 0) )) { 
      # If this is the first year a new cluster circulated, create new sets
      # Initialize new sets by copying the old sets
      YYs = rbind(YYs, YYs)
      # Then add exposure to the new strain to each of the newly branched sets
      YYs[(n.past.sets+1):(2^challenge.strain), challenge.strain] = 1
      #apply(YYs, 1, function(x) which(x > 0))
      
      # Update the fraction of people with past exposure to the new cluster
      frac.prior.exposure[challenge.strain, ] = p.infected[yy,]*OAS.free.frac
    }else{
      # If the challenge strain is part of a cluster that has previously circulated  
      # Update the fraction of people with past exposure to the strain of interest
      frac.prior.exposure[challenge.strain, ] = 
        frac.prior.exposure[challenge.strain, ] + (1-frac.prior.exposure[challenge.strain, ])*p.infected[yy,]*OAS.free.frac
      # Add the fraction previously exposed to (1-p(already exposed))*p(new exposure)
    }
    
    # Account for aging
    frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
    frac.prior.exposure[,1] = 0
    
  }
  
  
  # __________________ Age-specific incidence ____________________ #
  demog.yrs = as.character(c(rep(1980, 1980-start.year+1), 1981:max(years)))
  denom = rowSums(p.infected*demog[demog.yrs, ]); denom[which(denom == 0)] = 1
  outs = p.infected*demog[demog.yrs, ]/denom
  rownames(outs) = rownames(p.infected)
  outs
}



## Test
# simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 1.2, vv.hat = .5, tau.hat = .7, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, fraction.vaccinated = fraction.vaccinated, Hm = .5, wm.in = wm.H1, wo.in = wo.H1)
# simulate.incidence(alpha.hat = .2, AA.hat = .1, R0.hat = 2, vv.hat = .5, tau.hat = .7, Hm = .5,  demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, fraction.vaccinated = fraction.vaccinated, wo.in = wo.H3, wm.in = wm.H3)




############################################################
#____________ Define likelihood function  ______________#
#
############################################################
nll.prof = function(pars, prof.par, prof.par.value, demog, clusters, dat.in, st.yr, end.yr, frac.vaccinated, wm, wo, plot = FALSE){
  
  alpha.hat = ifelse(prof.par == 'alpha.hat', prof.par.value, pars['alpha.hat'])
  AA.hat = ifelse(prof.par == 'AA.hat', prof.par.value, pars['AA.hat']) 
  R0.hat = ifelse(prof.par == 'R0.hat', prof.par.value, pars['R0.hat'])
  vv.hat = ifelse(prof.par == 'vv.hat', prof.par.value, pars['vv.hat'])
  tau.hat = ifelse(prof.par == 'tau.hat', prof.par.value, pars['tau.hat'])
  Hm = ifelse(prof.par == 'Hm', prof.par.value, pars['Hm'])
  
  if(any(c(alpha.hat, AA.hat, R0.hat, vv.hat, tau.hat, Hm) < 0) | Hm > 1) {
    nll.yearly = rep(100000, 10) # Strong penalty for any parameter values that go negative
  }else{
  
  # __________________ Use model to generate expected probabilities 
  #___________________ that any case falls in a given birth year ____________________ #  
  sim = simulate.incidence(alpha.hat, AA.hat, R0.hat, vv.hat, tau.hat, Hm, demog = demog, clusters = clusters, start.year = st.yr, end.year = end.yr, fraction.vaccinated = frac.vaccinated, wo.in = wo, wm.in = wm)
  
  
  #___________________ format imprinting weights ____________________ # 
  fit.years = rownames(dat.in)
  # wm.in = wo.in = dat.in*0 # Initalize
  # for(ii in 1:nrow(dat.in)){
  #   yr = as.numeric(fit.years[ii])
  #   wm.in[ii, ] = as.numeric(wm[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
  #   wo.in[ii, ] = as.numeric(wo[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
  # }
  
  
  #___________________ calculate neg log lk ____________________ #  
  fit.years = rownames(dat.in)
  p.hat = sim[fit.years, ]
  nll.yearly = vector('numeric', length(fit.years))
  for(ii in 1:length(fit.years)){
    if( sum(p.hat[ii, ]) == 0 ){
      nll.yearly[ii] = 500  # Penalize if no outbreak occurs in the year of interest
    }else{
      nll.yearly[ii] = -dmultinom(x = as.numeric(dat.in[ii,]), prob = as.numeric(p.hat[ii,]), log = TRUE)
    }
  }
  }
  
  sum(nll.yearly)
}



## Test
load('Profs_H1_Ht.RData')
#nll.prof(pars = c(H1.OAS.Ht.profs[[6]]$par), prof.par = 'Hm', prof.par.value = .9583333, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, frac.vaccinated = fraction.vaccinated, wm = wm.H1, wo = wo.H1, plot = TRUE)
# nll.prof(pars = c(H1.OAS.Ht.profs[[6]]$par), prof.par = 'Hm', prof.par.value = Hm.vals.H1[6], demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, frac.vaccinated = fraction.vaccinated, wm = wm.H1, wo = wo.H1, plot = TRUE)
# 
# load('H1ests_vv_OAS_Ht2017-08-02.RData')
# nll.prof(pars = H1.ests.vv.OAS.Ht$par, prof.par = 'Hm', prof.par.value = .99998, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, frac.vaccinated = fraction.vaccinated, wm = wm.H1, wo = wo.H1, plot = TRUE)


nll.prof(pars = H1.OAS.Ht.profs[[6]]$par, prof.par = 'Hm', prof.par.value = Hm.vals.H1[6], demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, frac.vaccinated = fraction.vaccinated, wm = wm.H1, wo = wo.H1, plot = TRUE)


# # # 
# # # 
# 
# 




############################################################
#_________________ Calculate Hm profiles  ___________________#

#
print('Start H1N1 prof OAS Imprinting Transmission')
Sys.time()
#H1N1

# start a parallel cluster
library(parallel)
no_cores <- detectCores() - 3
cl <- makeCluster(no_cores)
# ## export relevant variables
clusterExport(cl, c('epi_final_size', 'simulate.incidence', 'nll.prof', 'demog', 'H1.clusters', 'H3.clusters', 'H1.inputs', 'H3.inputs', 'fraction.vaccinated', 'wm.H1', 'wm.H3', 'wo.H1', 'wo.H3', 'waifw.disaggregate', 'get.p.tns'))
# 
# Hm.vals.H1 = seq(.85, 1, length = 4*no_cores)
# 
# # calc profs
# H1.OAS.Ht.profs = parLapply(cl = cl, Hm.vals.H1,
#                             fun = function(pp.val) optim(par = c(alpha.hat = .5, AA.hat = 2.88, R0.hat = 3.3,
#                             vv.hat = .4, tau.hat = .6), fn = nll.prof, prof.par = 'Hm',
#                             prof.par.value = pp.val,  demog = demog, clusters = H1.clusters,
#                             dat.in = H1.inputs, st.yr = 1977, end.yr = 2014,
#                             frac.vaccinated = fraction.vaccinated, wm = wm.H1,
#                             wo = wo.H1, control = list(maxit = 1000)))
# 
# save(H1.OAS.Ht.profs, Hm.vals.H1, file = 'Profs_H1_Ht.RData')
# H1.OAS.HT.prof.vals = sapply(H1.OAS.Ht.profs, FUN = function(xx) xx$value)
# save(H1.OAS.Ht.profs, H1.OAS.HT.prof.vals, Hm.vals.H1, file = paste('Profs_H1_Ht', Sys.Date(), '.RData', sep = ''))
load('Profs_H1_Ht2017-08-04.RData')

Sys.time()
print('Start H3N2 prof valculation OAS Ht')
Hm.vals.H3 = seq(.6, 1, length = 2*no_cores)

H3.OAS.Ht.profs = parLapply(cl = cl, Hm.vals.H3,
                            fun = function(pp.val) optim(par = c(alpha.hat = .4, AA.hat = 3.27, R0.hat = 3.71,
                                                                 vv.hat = .23, tau.hat = .57), fn = nll.prof, prof.par = 'Hm',
                                                         prof.par.value = pp.val,  demog = demog, clusters = H3.clusters,
                                                         dat.in = H3.inputs, st.yr = 1968, end.yr = 2014,
                                                         frac.vaccinated = fraction.vaccinated, wm = wm.H3,
                                                         wo = wo.H3, control = list(maxit = 1000)))
save(H3.OAS.Ht.profs, Hm.vals.H3, file = 'Profs_H3_Ht.RData')
H3.OAS.HT.prof.vals = sapply(H3.OAS.Ht.profs, FUN = function(xx) xx$value)
save(H3.OAS.Ht.profs, H3.OAS.HT.prof.vals, Hm.vals.H3, file = 'Profs_H3_Ht.RData')


stopCluster(cl)

print('OAS Ht profs done')




## Analyze
load('Profs_H1_Ht2017-08-04.RData')
load('H1ests_vv_OAS_Ht2017-08-02.RData')
Hm.vals.H1
H1.OAS.HT.prof.vals - H1.ests.vv.OAS.Ht$value
plot(Hm.vals.H1, H1.OAS.HT.prof.vals - H1.ests.vv.OAS.Ht$value)
abline(h = 2)

load('Profs_H3_Ht.RData')
load('H3ests_vv_OAS_Ht2017-08-03.RData')
Hm.vals.H3
H3.OAS.HT.prof.vals - H3.ests.vv.OAS.Ht$value
plot(Hm.vals.H3, H3.OAS.HT.prof.vals - H3.ests.vv.OAS.Ht$value, ylim = c(-50, 50))
abline(v = H3.ests.vv.OAS.Ht$par['Hm'], col = 'blue')
abline(h = 2)


## Look at profiles against redo of optimization:
## Analyze
load('Profs_H1_Ht2017-08-04.RData')
load('Ht_ests_redo.RData')
Hm.vals.H1
H1.OAS.HT.prof.vals - H1.ests.vv.OAS.Ht.redo$value
plot(Hm.vals.H1, H1.OAS.HT.prof.vals - H1.ests.vv.OAS.Ht.redo$value, xlab = 'Hm value', ylab = 'Neg log lk difference')
abline(h = c(0,2))

load('Profs_H3_Ht.RData')
load('H3ests_vv_OAS_Ht2017-08-03.RData')
Hm.vals.H3
H3.OAS.HT.prof.vals - H3.ests.vv.OAS.Ht.redo$value
plot(Hm.vals.H3, H3.OAS.HT.prof.vals - H3.ests.vv.OAS.Ht$value, xlab = 'Hm value', ylab = 'Neg log lk difference')
abline(v = H3.ests.vv.OAS.Ht$par['Hm'], col = 'blue')
#abline(v = Hm.vals.H3[16], col = 'blue')
abline(h = 2)
abline(h = 0)

# prof.vals = 0
# ## Recalculate nlls using the best fit profile pars
# for(ii in 1:length(H1.OAS.Ht.profs)){
#   prof.vals[ii] = nll.prof(pars = H1.OAS.Ht.profs[[ii]]$par, prof.par = 'Hm', prof.par.value = Hm.vals.H1[ii], demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, frac.vaccinated = fraction.vaccinated, wm = wm.H1, wo = wo.H1, plot = TRUE)
# }

## Try re-running the full optimization with the best fit prof values


