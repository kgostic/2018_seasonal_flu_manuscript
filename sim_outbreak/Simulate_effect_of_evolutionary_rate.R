rm(list = ls())
setwd('~/Dropbox/R/2018_seasonal_flu/')
source('epi_final_size.R')


#### Simulate the age-specific incidence of seasonal influenza
#### Set a fixed R0, and assume year-to-year antigenic advance at different rates
#### Track the proportion of each age group with past exposure to all clusters {1, 2, ... n} in the set Y
####   Assume cross-protection decays exponentially with antigenic distance between strains
#### Use model formulation from http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002741



#_________ SET INPUT PARAMETERS for H1N1  __________# 
# Source: http://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1002741.s006&type=supplementary
# R0.hat = vary from 1.2 to 2.2
# AA.hat = 0.34
# alpha.hat = 0.25
# nu.hat = 0.06
#####################################################




#_________ Define demographic vectors __________# 
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



#_________ Define antigenic distance between strains circulating in simulated years __________#
get.aa.input = function(years, rate){
  
  antigenic.advance = seq(0, 0+rate*(length(years)-1), by = rate)

  data.frame('year' = years, 
                      'jump.year' = rep(1, length(years)),
                         'cluster' = antigenic.advance)
}


#_________ Define a funciton to calculate the prob of transmission   __________#
#_________       transmission given past exposure to a given set     __________#

# sets is a binary vector indicating which clusters are contained in the set
# challenge.value gives the index of the challenge strain
get.p.tns = function(set, challenge.value, all.values, AA.hat, alpha.hat){
  if(any(set == 1)){ # If any prior exposures
    in.set = which(set == 1)
    distances = abs(challenge.value - all.values[in.set])
    max(0, min(1-AA.hat*exp(-alpha.hat*distances)) )
  }else{ # Else, no prior exposures and fully transmisssible
    1
  }
}















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






############################################################
#____________ Simulation function  ______________#
#
############################################################
simulate.incidence = function(alpha.hat, AA.hat, R0.hat, demog, clusters = H1.clusters, start.year = 1977, end.year = 2014){
  
  # Use epi_final_size function
  # Initialize at 1977 and itterate forward to 2014 (the last year where we have data)
  years = start.year:end.year
  p.infected = matrix(NA, nrow = length(years), ncol = 100, dimnames = list(years, 0:99)) # Declare output storage
  
  ## Initialize a matrix to keep track of what fraction of each age group has previous exposure to each strain
  frac.prior.exposure = matrix(0, nrow = length(clusters$cluster), ncol = 100); colnames(frac.prior.exposure) = paste('age.', 0:99, sep = ''); rownames(frac.prior.exposure) = paste('cluster.', clusters$cluster, sep = '')
  ## Rows define ages, columns define sets and entries the define the fraction of the population in age group a with previous exposure to strain X
  
  
  # These matrices track which strains belong in each set
  # Each row is a different possible set
  # Columns take values 0 or 1 to indicate whether the set contains a strain from each cluster
  YYs = matrix(0, nrow = 1, ncol = length(clusters$cluster), dimnames = list(NULL, paste('cluster', clusters$cluster)))
  # First row defines the null set
  # Subsequent rows will be added within the simulation
  
  
  #_________ Caculate final sizes __________#
  # 1977 Step
  # Calculate final size, assume everyone is susceptible initially
  p.infected[1, ] = epi_final_size(r0 = 3, contact_matrix = waifw.disaggregate[,,1], demography_vector = demog[1,], prop_suscep = rep(1, 100))
  # Optional plot of results
  #plot(0:99, p.infected[1,]*demog[1,], main = years[1], ylab = 'Desnity of cases', xlab = 'age', ylim = c(0, .02))
  
  # Update the sets of possible past exposures
  YYs = rbind(YYs, YYs)
  # Then add the new strain to the set of possible xposures
  YYs[2, 1] = 1
  # Update the fraction of people with past exposure to the strain of interest
  frac.prior.exposure[1, ] = p.infected[1,]
  # Account for aging
  frac.prior.exposure[,2:100] = frac.prior.exposure[,1:99]
  frac.prior.exposure[,1] = 0
  
  
  
  # 1978 - 2017 epidemics
  for(yy in 2:length(years)){
    past.clusters = unique(clusters$cluster[1:(yy-1)]) # Vector of all clusters that have circulated previously
    n.past.clusters = length(past.clusters)
    n.past.sets = 2^(n.past.clusters-1)
    challenge.strain = clusters$cluster[yy]
    
    #___________ Calculate the fraction susceptible _________#
    #S_Ys is a matrix that gives the fraction of each age group with past exposure to set Y
    S_Ys = matrix(NA, n.past.sets+1, 100)
    
    for(ss in 1:nrow(YYs)){
      ## Multiply across all probabilities of having been previously exposed to exactly each set
      ## I.E, to each strain in the set, and of having not been previously exposed to each strain not in the set
      #S_Ys[ss, ] = apply((YYs[ss, ]*frac.prior.exposure+(1-YYs[ss, ])*(1-frac.prior.exposure)), 2, prod)
      S_Ys[ss, ] = exp(colSums(log( YYs[ss, ]*frac.prior.exposure+(1-YYs[ss, ])*(1-frac.prior.exposure) )))
    }
    
    #if(any(round(S_Ys, 6) != round(S_Ys.2, 6))) warning('New calculation not working')
    #if(!any(round(S_Ys, 6) != round(S_Ys.2, 6))) print('New calculation works!')
    
    if(any(round(colSums(S_Ys, na.rm = TRUE), 10) != 1)) warning('Faulty S_Y calculation, fractions do not sum to 1')
    
    
    
    # Get the fraction of each age group that can transmit the challenge strain
    p.tns = apply(YYs, 1, function(ss) get.p.tns(ss, challenge.strain, all.values = clusters$cluster, AA = AA.hat, alpha = alpha.hat))
    frac.susceptible = as.vector(t(S_Ys) %*% p.tns) # Gives a vector of the probability that each age group will transmit, given exposure histories
    
    

      p.infected[yy, ] = epi_final_size(r0 = R0.hat, contact_matrix = waifw.disaggregate[,,yy], demography_vector = demog[yy,], prop_suscep = frac.susceptible)
     

      # If this is the first year a new cluster circulated, create new sets
      # Initialize new sets by copying the old sets
      YYs = rbind(YYs, YYs)
      # Then add exposure to the new strain to each of the newly branched sets
      YYs[(n.past.sets+1):(2^(n.past.sets)), n.past.clusters] = 1
      #apply(YYs, 1, function(x) which(x > 0))
      
      # Update the fraction of people with past exposure to the new cluster
      frac.prior.exposure[challenge.strain, ] = p.infected[yy,]
    
    
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
simulate.incidence(alpha.hat = .2, AA.hat = .7, R0.hat = 1.2, demog = demog, clusters = get.aa.input(years = 2000:2020, rate = .9), start.year = 2000, end.year = 2020)
#simulate.incidence(alpha.hat = .2, AA.hat = .1, R0.hat = 2, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014)




############################################################
#____________ Define likelihood function  ______________#
#
############################################################
nll = function(pars, demog, clusters, dat.in, st.yr, end.yr, plot = FALSE){
  
  alpha.hat = pars['alpha.hat']  
  AA.hat = pars['AA.hat']  
  R0.hat = pars['R0.hat']  
  
  # __________________ Use model to generate expected probabilities 
  #___________________ that any case falls in a given birth year ____________________ #  
  sim = simulate.incidence(alpha.hat, AA.hat, R0.hat, demog = demog, clusters = clusters, start.year = st.yr, end.year = end.yr)
  
  
  
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
    if(plot == TRUE){
      plot(0:99, p.hat[ii, ], xlab = 'Age', ylab = 'Fraction of cases', main = paste(fit.years[ii], alpha.hat, AA.hat, R0.hat), ylim = c(0, .05))
      points(0:99, dat.in[ii, ]/sum(dat.in[ii, ]), col = 'blue')
    }
  }
  
  sum(nll.yearly)
}

## Test
#nll(pars = c(alpha.hat = .2, AA.hat = .1, R0.hat = 2), demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, plot = TRUE)
#nll(pars = c(alpha.hat = .2, AA.hat = .1, R0.hat = 2), demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, plot = TRUE)




############################################################
#_________________ Estimate par values  ___________________#
#
# H1N1
H1.ests = optim(par = c(alpha.hat = .1, AA.hat = .6, R0.hat = 2), fn = nll, demog = demog, clusters = H1.clusters, dat.in = H1.inputs, st.yr = 1977, end.yr = 2014, lower = c(0, 0, 0), upper = c(10, 10, 10), method = 'L-BFGS-B', plot = TRUE)

beepr::beep(); H1.ests
#
#save(H1.ests, file = paste('H1ests_', Sys.Date(), ".RData", sep = '' ))
#load('H1ests_2017-07-27.RData')
#
#
#
# H3N2
# H3.ests = optim(par = c(alpha.hat = .1, AA.hat = .6, R0.hat = 2), fn = nll, demog = demog, clusters = H3.clusters, dat.in = H3.inputs, st.yr = 1968, end.yr = 2014, lower = c(0, 0, 0), upper = c(10, 10, 10), method = 'L-BFGS-B', plot = TRUE)
# # 
# beepr::beep(); H3.ests
# #
# save(H3.ests, file = paste('H3ests_', Sys.Date(), ".RData", sep = '' ))
#
############################################################



#####################################################################
#____ Use best par values to generate age-incidence curves  ________#
## H1N1
pdf('null_model.pdf', height = 12)
par(mfrow = c(4, 2))
H1.simulation = simulate.incidence(alpha.hat = H1.ests$par['alpha.hat'], AA.hat = H1.ests$par['AA.hat'], R0.hat = H1.ests$par['R0.hat'], demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014)
# Plot simulated results and data
dat.years = rownames(H1.inputs)
for(ii in 1:nrow(H1.simulation)){
  yr = rownames(H1.simulation)[ii]
  plot(0:99, H1.simulation[ii, ], xlab = 'Age', ylab = 'Fraction of cases', main = rownames(H1.simulation)[ii], ylim = c(0, .05))
  # For years with more than 100 observed data points, plot the data too
  if(yr %in% dat.years){
    if(sum(H1.inputs[yr, ]) > 100){
      points(0:99, H1.inputs[rownames(H1.simulation)[ii], ]/sum(H1.inputs[rownames(H1.simulation)[ii], ]), col = 'blue')
    }}
}
dev.off()



## H3N2
H3.simulation = simulate.incidence(alpha.hat = H3.ests$par['alpha.hat'], AA.hat = H3.ests$par['AA.hat'], R0.hat = H3.ests$par['R0.hat'], demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014)
# Plot simulated results and data
dat.years = rownames(H3.inputs)
for(ii in 1:nrow(H3.simulation)){
  yr = rownames(H3.simulation)[ii]
  plot(0:99, H3.simulation[ii, ], xlab = 'Age', ylab = 'Fraction of cases', main = rownames(H3.simulation)[ii], ylim = c(0, .05))
  # For years with more than 100 observed data points, plot the data too
  if(yr %in% dat.years){
    if(sum(H3.inputs[yr, ]) > 100){
      points(0:99, H3.inputs[rownames(H3.simulation)[ii], ]/sum(H3.inputs[rownames(H3.simulation)[ii], ]), col = 'blue')
    }}
}

dev.off()
