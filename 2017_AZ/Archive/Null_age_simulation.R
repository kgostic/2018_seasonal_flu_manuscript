rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Kucharski_SIRfinalsize/R/epi_final_size.R')

#### Simulate the age-specific incidence of seasonal influenza
#### Begin simulation in 1968 for H3N2 and in 1977 for H1N1
#### Track the proportion of each age group with past exposure to all clusters {1, 2, ... n} in the set Y
#### Assume all seasonal variants in the same cluster are antigenically identical, and decreasing cross immunity between strains in increasingly distant clusters.
#### Use model formulation from http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002741
#### Assume cluster jumps as in http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s005&type=supplementary

## !!!!!!!!!!! Marks things that need to be polished


#_________ SET INPUT PARAMETERS for H1N1  __________# 
# Source: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s006&type=supplementary
R0 = 1.2
AA = 0.34
alpha = 0.25
nu = 0.06
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



#_________ Define sets of strains to which prior exposure is possible __________#
# n.pos.sets.H1 = 2^c(0, H1.clusters$cluster[-nrow(H1.clusters)]); names(n.pos.clust.H1) = 1977:2017
# n.pos.sets.H3 = 2^c(0, H3.clusters$cluster[-nrow(H3.clusters)]); names(n.pos.clust.H3) = 1968:2017
# # Initialize a list of all possible sets
# YYs.H1 = vector('list', 1) # First set of possible past exposures is the null
# YYs.H3 = vector('list', 1)
# names(YYs.H1) = 1977
# names(YYs.H3) = 1968

# These matrices track which strains belong in each set
# Each row is a different possible set
# Columns take values 0 or 1 to indicate whether the set contains a strain from each cluster
YYs.H1 = matrix(0, nrow = 1, ncol = max(H1.clusters$cluster), dimnames = list(NULL, paste('cluster', 1:max(H1.clusters$cluster))))
# First row defines the null set
# Subsequent rows will be added within the simulation


## Define a funciton to calculate the prob of transmission given prior exposure
# sets is a binary vector indicating which clusters are contained in the set
# challenge.value gives the index of the challenge strain
get.p.tns = function(set, challenge.value, AA, alpha){
  if(any(set == 1)){ # If any prior exposures
  distances = abs(challenge.value - which(set == 1))
  max(0, min(1-AA*exp(-alpha*distances)))
  }else{ # Else, no prior exposures and fully transmisssible
  return(1)
  }
}


## Define a function to calculate Q(a), the fraction of each age group that should be able to transmit
# get.Qas = function(sets, challenge.value){
#   options(warn = -1)
#   set.indices = 1:length(sets)
#   sapply(0:99, function(ages) sum(sapply(set.indices, function(x) p.tns(set.index = x, sets = sets, challenge.value =1)*S_Y[ages+1, x]), na.rm = TRUE))
# }

## Initialize a matrix to keep track of what fraction of each age group has previous exposure to each strain
frac.prior.exposure = matrix(0, nrow = max(H1.clusters$cluster), ncol = 100); colnames(frac.prior.exposure) = paste('age.', 0:99, sep = ''); rownames(frac.prior.exposure) = paste('cluster.', 1:max(H1.clusters$cluster), sep = '')
## Rows define ages, columns define sets and entries the define the fraction of the population in age group a with previous exposure to strain X






####################################################
#       LIST OF VARIABLES IN SIMULATION
#
# p.infected.H1 - Matrix of probabilities that a person of any age from 0:99 
#                 became infected in the season of interest.
# YYs.H1 - Binary matrix of sets. Each row represents a possible set of past exposures.
#                 Each column represents a cluster, 0 or 1 indicates whether the cluster
#                 belongs to the set.
# frac.prior exposure - Matrix of probabilities that a person of a given age has past 
#                 exposure to each individual cluster. Updated each itteration.
#
#
#
#







##############################################
# __________ BEGIN H1 SIMULATION   __________#
# Use epi_final_size function
# Initialize at 1977 and itterate forward to 2017
years.H1 = 1977:2017
p.infected.H1 = matrix(NA, nrow = length(years.H1), ncol = 100, dimnames = list(years.H1, 0:99)) # Declare output storage


#_________ Caculate final sizes __________#
# 1977 Step
# Calculate final size, assume everyone is susceptible initially
p.infected.H1[1, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = rep(1, 100))
# Optional plot of results
plot(0:99, p.infected.H1[1,]*demog[1,], main = years.H1[1], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))

# Update the sets of possible past exposures
YYs.H1 = rbind(YYs.H1, YYs.H1)
# Then add exposure to the new strain to each of the new sets
YYs.H1[2, 1] = 1
# Update the fraction of people with past exposure to the strain of interest
frac.prior.exposure[1, ] = p.infected.H1[1,]
# Account for aging
frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
frac.prior.exposure[,1] = 0



# 1978 - 2017 epidemics
for(yy in 2:length(years.H1)){
  past.clusters = unique(H1.clusters$cluster[1:yy-1]) # Vector of all clusters that have circulated previously
  n.past.clusters = length(past.clusters)
  n.past.sets = 2^n.past.clusters
  challenge.strain = H1.clusters$cluster[yy]

#___________ Calculate the fraction susceptible _________#
#S_Ys is a matrix that gives the fraction of each age group with past exposure to set Y
S_Ys = matrix(NA, n.past.sets, 100)
  
for(ss in 1:nrow(YYs.H1)){
  ## Multiply across all probabilities of having been previously exposed to each strain in the set, and of having not been previously exposed to each strain not in the set
  #S_Ys[ss, ] = apply((YYs.H1[ss, ]*frac.prior.exposure+(1-YYs.H1[ss, ])*(1-frac.prior.exposure)), 2, prod)
  S_Ys[ss, ] = exp(colSums(log( YYs.H1[ss, ]*frac.prior.exposure+(1-YYs.H1[ss, ])*(1-frac.prior.exposure) )))
}

#if(any(round(S_Ys, 6) != round(S_Ys.2, 6))) warning('New calculation not working')
#if(!any(round(S_Ys, 6) != round(S_Ys.2, 6))) print('New calculation works!')

if(any(round(colSums(S_Ys, na.rm = TRUE), 10) != 1)) warning('Faulty S_Y calculation, fractions do not sum to 1')



# Get the fraction of each age group that can transmit the challenge strain
p.tns = apply(YYs.H1, 1, function(ss) get.p.tns(ss, challenge.strain, AA = AA, alpha = alpha))
frac.susceptible = as.vector(t(S_Ys) %*% p.tns) # Gives a vector of the probability that each age group will transmit, given exposure histories
  

  if(years.H1[yy] <= 1980){ # 1980 is the first year in which demographic information was available, so use demog[1, ] vector for all years from 1977:1980
    p.infected.H1[yy, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = frac.susceptible)
    plot(0:99, p.infected.H1[yy,]*demog[yy,], main = years.H1[yy], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
    
  }else{
    rr = 2
    p.infected.H1[yy, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,rr], demography_vector = demog[rr,], prop_suscep = frac.susceptible)
    plot(0:99, p.infected.H1[yy,]*demog[yy,], main = years.H1[yy], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
    rr = rr+1
  }

# Update the sets of possible past exposures
if(challenge.strain > max( which(colSums(YYs.H1) > 0) )) { 
  # If this is the first year a new cluster circulated, create new sets
  # Initialize new sets by copying the old sets
  YYs.H1 = rbind(YYs.H1, YYs.H1)
  # Then add exposure to the new strain to each of the newly branched sets
  YYs.H1[(n.past.sets+1):(2^challenge.strain), challenge.strain] = 1
  #apply(YYs.H1, 1, function(x) which(x > 0))
  
  # Update the fraction of people with past exposure to the new cluster
  frac.prior.exposure[challenge.strain, ] = p.infected.H1[yy,]
}else{
  # If the challenge strain is part of a cluster that has previously circulated  
  # Update the fraction of people with past exposure to the strain of interest
  frac.prior.exposure[challenge.strain, ] = 
    frac.prior.exposure[challenge.strain, ] + (1-frac.prior.exposure[challenge.strain, ])*p.infected.H1[yy,]
  # Add the fraction previously exposed to (1-p(already exposed))*p(new exposure)
}

# Account for aging
frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
frac.prior.exposure[,1] = 0
 
}


cols = rainbow(50)[-(42:50)]
plot(0:99, p.infected.H1[1, ]/sum(p.infected.H1[1,]), col = cols[1], type = 'l', xlab = 'age', ylab = 'Fraction of infections', ylim = c(0, .02), main = 'H1N1 predictions')
for(ii in 2:41){
lines(0:99, p.infected.H1[ii, ]/sum(p.infected.H1[ii,]), col = cols[ii], lwd = 4*seq(1, 1/4, length = 41)[ii])
}


barplot(frac.prior.exposure, col = rainbow(11))











#_________ SET INPUT PARAMETERS for H3N2  __________# 
# Source: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s006&type=supplementary
R0 = 1.63
AA = 1.2
alpha = 0.17
nu = 0.24
#####################################################

## Initialize a matrix to keep track of what fraction of each age group has previous exposure to each cluster
frac.prior.exposure = matrix(0, nrow = max(H3.clusters$cluster), ncol = 100); colnames(frac.prior.exposure) = paste('age.', 0:99, sep = ''); rownames(frac.prior.exposure) = paste('cluster.', 1:max(H3.clusters$cluster), sep = '')
## Rows define ages, columns define sets and entries the define the fraction of the population in age group a with previous exposure to strain X
YYs.H3 = matrix(0, nrow = 1, ncol = max(H3.clusters$cluster), dimnames = list(NULL, paste('cluster', 1:max(H3.clusters$cluster))))


##############################################
# __________ BEGIN H3 SIMULATION   __________#
years.H3 = 1968:2017
p.infected.H3 = matrix(NA, nrow = length(years.H3), ncol = 100, dimnames = list(years.H3, 0:99)) # Declare output storage


#_________ Caculate final sizes __________#
# Use epi_final_size function
# Initialize at 1968 and itterate forward to 2017


# 1968 Step
# Calculate final size, assume everyone is susceptible initially
p.infected.H3[1, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = rep(1, 100))
# Optional plot of results
plot(0:99, p.infected.H3[1,]*demog[1,], main = years.H3[1], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))

# Update the sets of possible past exposures
YYs.H3 = rbind(YYs.H3, YYs.H3)
# Then add exposure to the new strain to each of the new sets
YYs.H3[2, 1] = 1
# Update the fraction of people with past exposure to the strain of interest
frac.prior.exposure[1, ] = p.infected.H3[1,]
# Account for aging
frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
frac.prior.exposure[,1] = 0






# 1978 - 2017 epidemics
for(yy in 2:length(years.H3)){
  past.clusters = unique(H3.clusters$cluster[1:yy-1]) # Vector of all clusters that have circulated previously
  n.past.clusters = length(past.clusters)
  n.past.sets = 2^n.past.clusters
  challenge.strain = H3.clusters$cluster[yy]
  
  #___________ Calculate the fraction susceptible _________#
  #S_Ys is a matrix that gives the fraction of each age group with past exposure to set Y
  S_Ys = matrix(NA, n.past.sets, 100)
  
  
  
  
  
  
  for(ss in 1:nrow(YYs.H3)){
    ## Multiply across all probabilities of having been previously exposed to each strain in the set, and of having not been previously exposed to each strain not in the set
    #S_Ys[ss, ] = apply((YYs.H3[ss, ]*frac.prior.exposure+(1-YYs.H3[ss, ])*(1-frac.prior.exposure)), 2, prod)
    S_Ys[ss, ] = exp(colSums(log( YYs.H3[ss, ]*frac.prior.exposure+(1-YYs.H3[ss, ])*(1-frac.prior.exposure) )))
  }
  
  if(any(round(colSums(S_Ys, na.rm = TRUE), 10) != 1)) warning('Faulty S_Y calculation, fractions sum to greater than 1')

  
  
  # Get the fraction of each age group that can transmit the challenge strain
  p.tns = apply(YYs.H3, 1, function(ss) get.p.tns(ss, challenge.strain, AA = AA, alpha = alpha))
  frac.susceptible = as.vector(t(S_Ys) %*% p.tns) # Gives a vector of the probability that each age group will transmit, given exposure histories
  
  
  if(years.H3[yy] <= 1980){ # 1980 is the first year in which demographic information was available, so use demog[1, ] vector for all years from 1977:1980
    p.infected.H3[yy, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = frac.susceptible)
    plot(0:99, p.infected.H3[yy,]*demog[yy,], main = years.H3[yy], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
    
  }else{
    rr = 2
    p.infected.H3[yy, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,rr], demography_vector = demog[rr,], prop_suscep = frac.susceptible)
    plot(0:99, p.infected.H3[yy,]*demog[yy,], main = years.H3[yy], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
    rr = rr+1
  }
  
  # Update the sets of possible past exposures
  if(challenge.strain > max( which(colSums(YYs.H3) > 0) )) { 
    # If this is the first year a new cluster circulated, create new sets
    # Initialize new sets by copying the old sets
    YYs.H3 = rbind(YYs.H3, YYs.H3)
    # Then add exposure to the new strain to each of the new sets
    YYs.H3[(n.past.sets+1):(2^challenge.strain), challenge.strain] = 1
    #apply(YYs.H3, 1, function(x) which(x > 0))
    
    # Update the fraction of people with past exposure to the strain of interest
    frac.prior.exposure[challenge.strain, ] = p.infected.H3[yy,]
  }else{
    # If the challenge strain is part of a cluster that has previously circulated  
    # Update the fraction of people with past exposure to the strain of interest
    frac.prior.exposure[challenge.strain, ] = 
      frac.prior.exposure[challenge.strain, ] + (1-frac.prior.exposure[challenge.strain, ])*p.infected.H3[yy,]
    # Add the fraction previously exposed to (1-previously.exposed)*p(new exposure)
  }
  
  # Account for aging
  frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
  frac.prior.exposure[,1] = 0

  
}




cols = rainbow(60)[-(51:60)]
plot(0:99, p.infected.H3[1, ]/sum(p.infected.H3[1,]), col = cols[1], type = 'l', xlab = 'age', ylab = 'Fraction of infections', ylim = c(0, .1), main = 'H3N2 predictions')
for(ii in 2:50){
  lines(0:99, p.infected.H3[ii, ]/sum(p.infected.H3[ii,]), col = cols[ii], lwd = 4*seq(1, 1/4, length = 41)[ii])
}




save(demog, p.infected.H1, p.infected.H3, file = 'Age_outputs.RData')
