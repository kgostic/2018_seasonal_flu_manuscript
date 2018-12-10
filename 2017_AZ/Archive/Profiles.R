rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Kucharski_SIRfinalsize/R/epi_final_size.R')
source('Source_functions.R')

#### Simulate the age-specific incidence of seasonal influenza
#### Begin simulation in 1968 for H3N2 and in 1977 for H1N1
#### Track the proportion of each age group with past exposure to all clusters {1, 2, ... n} in the set Y
#### Assume all seasonal variants in the same cluster are antigenically identical, and decreasing cross immunity between strains in increasingly distant clusters.
#### Use model formulation from http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002741
#### Assume cluster jumps as in http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s005&type=supplementary




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


### Disaggregate the WAIFW matrix from 5-year age groups to single-year age groups
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
#  * Units are average number of contacts per day made by a person of age i with a person of age j
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

rowSums(H1.inputs) # Simulate up to 2014 when fitting the model

#____________________ Load data on USA vaccination ______________________#
raw.frac.vaccinated = read.csv('USA_vaccination_coverage.csv', skip = 1, header = TRUE)
fraction.vaccinated = matrix( c(raw.frac.vaccinated[,2], # age 0
                                rep(raw.frac.vaccinated[,3], 17), # age 1-17
                                rep(raw.frac.vaccinated[,4], 32), # age 18-49 
                                rep(raw.frac.vaccinated[,5], 15), # age 50-64
                                rep(raw.frac.vaccinated[,6], 35)),# age 65-99
                              nrow = length(1968:2017), ncol = 100, dimnames = list(1968:2017, 0:99))/100



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










#_________ SET INPUT PARAMETERS for H1N1  __________# 
#Source: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s006&type=supplementary

# load('H3ests_2017-08-15.RData')
# 
# R0.hat = 2.5 #H3.ests$par['R0.hat']
# AA.hat = .7 #H3.ests$par['AA.hat']
# alpha.hat = 1.7 #H3.ests$par['alpha.hat']
# vv.hat = .5
# tau.hat = .7
# clusters = H3.clusters
# start.year = 1968
# end.year = 2014
# WAIFW = waifw.disaggregate
# vaccinate = TRUE
# Fraction.vaccinated = fraction.vaccinated
# OAS = TRUE
# imprinting = FALSE
# wm.in = wm.H3
# wo.in = wo.H3

#####################################################





############################################################
#____________ Simulation function  ______________#
#
############################################################
simulate.incidence = function(alpha.hat, AA.hat, R0.hat, vv.hat, tau.hat, Hm = 1, demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = NULL){

  # Use epi_final_size function
  # Initialize at 1977 and itterate forward to 2014 (the last year where we have data)
  years = start.year:end.year
  wm.use = wm.in[as.character(years), ]
  wo.use = 1-wm.use
  total.n.clusters = max(clusters$cluster)
  frac.vaccinated = Fraction.vaccinated[as.character(years), ]
  p.infected = matrix(NA, nrow = length(years), ncol = 100, dimnames = list(years, 0:99)) # Declare output storage
  
  ## Initialize a matrix to keep track of what fraction of each age group has previous exposure to each strain
  frac.prior.exposure = matrix(0, nrow = max(clusters$cluster), ncol = 100); colnames(frac.prior.exposure) = paste('age.', 0:99, sep = ''); rownames(frac.prior.exposure) = paste('cluster.', 1:max(clusters$cluster), sep = '')
  ## Rows define ages, columns define sets and entries the define the fraction of the population in age group a with previous exposure to strain X
  
  
  #_________ Caculate final sizes __________#
  # 1977 Step
  # Calculate final size, assume everyone is susceptible initially
  
  # If considering vaccination, reduce the initial susceptible fraction by vaccine coverage*effectiveness
  if(vaccinate == TRUE & imprinting != TRUE){ 
  frac.transmissible = rep(1, 100) * (1 - frac.vaccinated[1,]*vv.hat)
  # Else if no vaccination, assume everyone is initially susceptible
  }else if (vaccinate != TRUE & imprinting == TRUE) {
  frac.transmissible = rep(1, 100) * (Hm*wm.use[1, ]+wo.use[1, ])
  }else if (vaccinate == TRUE & imprinting == TRUE){
  frac.transmissible = rep(1, 100) * (1 - frac.vaccinated[1,]*vv.hat) * (Hm*wm.use[1, ] + wo.use[1, ])
  }else{
  frac.transmissible = rep(1, 100)
  }
  #print(frac.transmissible)
  
  p.infected[1, ] = epi_final_size(r0 = R0.hat, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = frac.transmissible)
 # print(p.infected[1,])

  # Update the fraction of people with past exposure to the strain of interest
  frac.prior.exposure[1, ] = p.infected[1, ]
  # Account for aging
  frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
  frac.prior.exposure[,1] = 0
  
  # Assume no OAS in the first year, as no one has imprinted yet
  
  
  
  
  # 1978 - 2017 epidemics
  for(yy in 2:length(years)){
    past.clusters = unique(clusters$cluster[1:yy-1]) # Vector of all clusters that have circulated previously
    challenge.cluster = clusters$cluster[yy]
    
    #----------------------------------------------------------------------------------------------#
    #___________ Calculate the fraction of each age group with last exposure to cluster x _________#
    #______________________________________________________________________________________________#
    P_x = apply(frac.prior.exposure, 2, get_last_exposure_probs)
    # P_x is a matrix representing the probability that in individual of age i's last flu exposure was to cluster x
    #   Each column represents a different age
    #   Each row represents a different cluster to which the last exposure took place
    #        Frist row represents no prior exposure (null set)
    #        Second row represents probabilty that the last exposure was to cluster 1
    #        Third row represents probability that he last exposure was to clsuter 2
    #        etc.

    
    #----------------------------------------------------------------------------------------------#
    #___________ Calculate the fraction of each age group that can transmit the challenge _________#
    #___________                 strain, given their most recent exposure                 _________#
    #______________________________________________________________________________________________#
    p.tns = c(1, 1-AA.hat*exp(-alpha.hat* (1:total.n.clusters - challenge.cluster)) )
    if(any(p.tns < 0)){ p.tns[which(p.tns<0)] = 0} # minimum prob of 0
    # Vector giving the probability of transmission, assuming that the most recent exposure was to cluster x
    # First entry is always 1, as individuals with no previous flu exposures are assumed fully susceptible
    # Second entry represents probability of transmitting if the most recent exposure was to cluster 1
    # etc. according to conventions for rows of P_x
    
    #--------------------------------------------------------------------------------------------------#
    #___________ Calculate the overall probability that an individual of age i will transmit  _________#
    #___________         given that age group's exposure histories                            _________#
    #__________________________________________________________________________________________________# 
    if(vaccinate == TRUE & imprinting != TRUE){ # If the model considers vaccination
      frac.transmissible = as.vector(t(P_x) %*% p.tns) * (1 - frac.vaccinated[yy,] * vv.hat)
    }else if (vaccinate != TRUE & imprinting == TRUE) { # If imprinting but not vaccination is included
      frac.transmissible = as.vector(t(P_x) %*% p.tns) * (Hm*wm.use[yy, ]+wo.use[yy, ])
    }else if (vaccinate == TRUE & imprinting == TRUE) { # If both imprinting and vaccination included
      frac.transmissible = as.vector(t(P_x) %*% p.tns) * (1 - frac.vaccinated[yy,] * vv.hat) * (Hm*wm.use[yy, ]+wo.use[yy, ])
    }else{ # If neither vaccination nor imprinting included
    frac.transmissible = as.vector(t(P_x) %*% p.tns) 
    }
    
    # Vector of probabilities that each age group will transmit. Length 100 with entries representing ages 0:99
    

    #--------------------------------------------------------------------------------------------------#
    #___________            Approximate final outbreak size  in each age group                _________#
    #__________________________________________________________________________________________________# 
    rr = max(1, years[yy] - 1979) # Use the first row of the demography vector and first slice of WAIFW matrix until 1980
    p.infected[yy, ] = epi_final_size(r0 = R0.hat, contact_matrix = WAIFW[,,rr], demography_vector = demog[rr,], prop_suscep = frac.transmissible)
    
    
    
    
    # Update the sets of possible past exposures
    if(as.logical(clusters$jump.year[yy])) { 
      # If this is the first year a new cluster circulated, update the fraction of each age group that has now been exposed to the new cluster
      if(OAS == TRUE & any(p.tns[1:challenge.cluster] < tau.hat)){
        # If model assumes OAS, any any exposure history would prevent the formation of new antibodies against the challenge strain, then reduce the fraction of individuals who will add to their antibody repitoire appropriately
        
        OAS.clusters = which(p.tns[1:challenge.cluster] < tau.hat)
        first.exposure.probs = apply(frac.prior.exposure, 2, get_first_exposure_probs)
        if(length(OAS.clusters) > 1){
          OAS.frac = colSums(first.exposure.probs[OAS.clusters, ])
        }else{
          OAS.frac = first.exposure.probs[OAS.clusters, ]
        }
        frac.prior.exposure[challenge.cluster, ] = p.infected[yy,]*(1-OAS.frac)
      }else{
        
      # IF we assume no OAS  
      frac.prior.exposure[challenge.cluster, ] = p.infected[yy,]
      }
      

      
    }else{
    # If the challenge strain is part of a cluster that has previously circulated  
    # Update the fraction of people with past exposure to the strain of interest
      
      if(OAS == TRUE & any(p.tns[1:challenge.cluster] < tau.hat)){
        # If model assumes OAS, any any exposure history would prevent the formation of new antibodies against the challenge strain, then reduce the fraction of individuals who will add to their antibody repitoire appropriately
        
        OAS.clusters = which(p.tns[1:challenge.cluster] < tau.hat)
        first.exposure.probs = apply(frac.prior.exposure, 2, get_first_exposure_probs)
        
        if(length(OAS.clusters) > 1){
          OAS.frac = colSums(first.exposure.probs[OAS.clusters, ])
        }else{
          OAS.frac = first.exposure.probs[OAS.clusters, ] 
          # OAS frac gives the fraction of each age group that is prevented from generating new immune memory
          # 1-OAS frac gives the fraction of each age group that will generate new immune memory, given exposure to the new cluster
        }
        
      frac.prior.exposure[challenge.cluster, ] = 
        frac.prior.exposure[challenge.cluster, ] + (1-frac.prior.exposure[challenge.cluster, ])*p.infected[yy,]*(1-OAS.frac)
      # Add the fraction previously exposed to (1-p(already exposed))*p(new exposure)*p(OAS does not interfere with new immune memory)
      
      
      }else{ # If the model assumes no OAS
      
      frac.prior.exposure[challenge.cluster, ] = 
        frac.prior.exposure[challenge.cluster, ] + (1-frac.prior.exposure[challenge.cluster, ])*p.infected[yy,]
      }
    }
    
    # Account for aging
    frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
    frac.prior.exposure[,1] = 0
    
    # Optional plot
    if(plot == TRUE){ plot(0:99, p.infected[yy, ], main = years[yy])}
  }
  
  
  # __________________ Age-specific incidence ____________________ #
  demog.yrs = as.character(c(rep(1980, 1980-start.year+1), 1981:max(years)))
  inf.age.dist = p.infected*demog[demog.yrs, ] # Gives the proportion of the population that is age a, and was infected
  denom = rowSums(p.infected*demog[demog.yrs, ]); denom[which(denom == 0)] = 1 # Gives the total proportion of the population infected
  outs = inf.age.dist/denom # Gives the proportion of all cases that occurred in each age group
  rownames(outs) = rownames(p.infected)
  outs
}



## Test
# no.vaccine = simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, Hm = .5, tau.hat = .7, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = FALSE, imprinting = TRUE, wm.in = wm.H1)
# vaccin.not.effective = simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 0, tau.hat = .7, Hm = .5, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = TRUE, imprinting = TRUE, wm.in = wm.H1)
# vaccine = simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = TRUE)
  
# no.OAS =        simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, tau.hat = .5, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = FALSE, OAS = FALSE)
# OAS.no.radius = simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, tau.hat = 0, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = FALSE, OAS = TRUE)
# OAS =           simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, tau.hat = .7, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = TRUE, OAS = TRUE)
# 
# no.imprinting  =           simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, tau.hat = .7, Hm = .5, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = FALSE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H1)
# no.imprinting.effect  =    simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, tau.hat = .7, Hm = 1, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = FALSE, OAS = TRUE, imprinting = TRUE, wm.in = wm.H1)
# imprinting  =           simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 2.5, vv.hat = 1, tau.hat = .7, Hm = 0, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, plot = FALSE, Fraction.vaccinated = fraction.vaccinated, vaccinate = TRUE, OAS = TRUE, imprinting = TRUE, wm.in = wm.H1)
# 
# par(mfrow = c(2,1))
# plot(0:99, no.imprinting[1, ])
# lines(0:99, imprinting[1, ], col = 'red')
# barplot(wm.H1[rownames(H1.inputs), ])









############################################################
#____________ Define likelihood function  ______________#
#
############################################################
nll_prof = function(pars, Hm, demog, clusters, dat.in, st.yr, end.yr, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = NULL, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = FALSE){
  
  alpha.hat = pars['alpha.hat']  
  AA.hat = pars['AA.hat']  
  R0.hat = pars['R0.hat']
  if(vaccinate.switch == TRUE){
    vv.hat = pars['vv.hat'] 
  }else{
    vv.hat = NULL
  }
  if(OAS.switch == TRUE){
    tau.hat = pars['tau.hat'] 
  }else{
    tau.hat = NULL
  }

  
  if(any(c(alpha.hat, AA.hat, R0.hat, vv.hat, Hm, tau.hat) < 0)) {
    nll.yearly = rep(100000, 10) # Strong penalty for any parameter values that go negative
  }else{

  # __________________ Use model to generate expected probabilities 
  #___________________ that any case falls in a given birth year ____________________ #  
  sim = simulate.incidence(alpha.hat, AA.hat, R0.hat, vv.hat, tau.hat, Hm, demog = demog, clusters = clusters, start.year = st.yr, end.year = end.yr, WAIFW = WAIFW.in, vaccinate = vaccinate.switch, Fraction.vaccinated = vax.matrix, OAS = OAS.switch, imprinting = imprinting.switch, wm.in = wm.in)
  
  
  
  #___________________ calculate neg log lk ____________________ #  
  fit.years = rownames(dat.in)
  p.hat = sim[fit.years, ]
  nll.yearly = vector('numeric', length(fit.years))
  for(ii in 1:length(fit.years)){
    if( sum(p.hat[ii, ]) == 0 ){
      #print(ii)
      nll.yearly[ii] = 500  # Penalize if no outbreak occurs in the year of interest
    }else{
      nll.yearly[ii] = -dmultinom(x = as.numeric(dat.in[ii,]), prob = as.numeric(p.hat[ii,]), log = TRUE)
    }
    
    
    if(plot == TRUE ){
      cols = rainbow(length(fit.years))
      if(ii == 1){
    plot(0:99, p.hat[ii, ], xlab = 'Age', ylab = 'Fraction of cases', main = paste(round(R0.hat, 3), round(AA.hat, 3), round(alpha.hat, 3)), ylim = c(0, .05), col = cols[ii], type = 'l')
        points(0:99, dat.in[ii, ]/sum(dat.in[ii, ]), col = cols[ii])
      }else{
        lines(0:99, p.hat[ii, ], col = cols[ii])
     points(0:99, dat.in[ii, ]/sum(dat.in[ii, ]), col = cols[ii])
      }
      if(ii == length(fit.years))      text(60, .04, paste(sum(nll.yearly)))
    }
    
    
  }
  }
  
  sum(nll.yearly)
}


#nll_prof(pars = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7), Hm = .5, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H1)

#optim(par = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7), nll_prof, Hm = .5, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H1)
# 
# nll(pars = c(AA.hat = .7, alpha.hat = .6, R0.hat = 2.5, vv.hat = .5, tau.hat = .7, Hm = .5), demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = TRUE, wm.in = wm.H1)


