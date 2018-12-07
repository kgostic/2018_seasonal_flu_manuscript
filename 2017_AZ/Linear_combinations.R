## Linear combinations models discussed with Jamie

rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')
source('Kucharski_SIRfinalsize/R/epi_final_size.R')
source('Source_functions.R')
source('Baseline.R')

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
  valid = which(rowSums(out.mat) > 300)
  out.mat[valid, ]
}

## Reformatted data inputs
H1.inputs = reformat.data(H1_AZ)
H3.inputs = reformat.data(H3_AZ)

rowSums(H1.inputs) # Simulate up to 2014 when fitting the model
rowSums(H3.inputs)

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
# wm.H1 = wo.H1 = w1.formatted =  fraction.vaccinated*0 # Initalize
# for(ii in 1:nrow(fraction.vaccinated)){
#   yr = as.numeric(rownames(fraction.vaccinated))[ii]
#   wm.H1[ii, ] = as.numeric((w1+w2)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
#   wo.H1[ii, ] = as.numeric((w3+wn)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
#   w1.formatted[ii, ] = as.numeric((w1)[paste(yr, 'USA', sep = ''), as.character(yr:(yr-99))])
# }
# 
# save(w1.formatted, file = 'Group1_weights.RData')

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

# Use Adam's parameters for H1
H1.baseline = simulate.incidence(alpha.hat = .25, AA.hat = .34, R0.hat = 1.3, vv.hat = .06, tau.hat = 0.93, Hm = 1, demog = demog, clusters = H1.clusters, start.year = 1977, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = (waifw.disaggregate), plot = FALSE, vaccinate = FALSE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H1)

cols = rainbow(30)
plot(0:99, H1.baseline[1,], col = cols[1], type = 'l', ylim = c(0, .07))
for(ii in 2:nrow(H1.baseline)){lines(0:99, H1.baseline[ii, ], col = cols[ii])}


# Use Adam's parameters for H3
H3.baseline = simulate.incidence(alpha.hat = .17, AA.hat = 1.20, R0.hat = 1.7, vv.hat = .24, tau.hat = .93, Hm = 1, demog = demog, clusters = H3.clusters, start.year = 1968, end.year = 2014, Fraction.vaccinated = fraction.vaccinated, WAIFW = waifw.disaggregate, plot = FALSE, vaccinate = TRUE, OAS = TRUE, imprinting = FALSE, wm.in = wm.H3)

cols = rainbow(45)
plot(0:99, H3.baseline[1,], col = cols[1], type = 'l', ylim = c(0, .06))
for(ii in 2:nrow(H3.baseline)){lines(0:99, H3.baseline[ii, ], col = cols[ii])}









nll = function(pars, demog, clusters, dat.in, st.yr, end.yr, plot = FALSE, WAIFW.in = waifw.disaggregate, vax.matrix = fraction.vaccinated, wm.in = NULL, vaccinate.switch = TRUE, OAS.switch = TRUE, imprinting.switch = FALSE, mass.action = FALSE){
  
  if(st.yr == 1977){
    penalty = 800
  }else if (st.yr == 1968){
    penalty = 1300
  }else{
    penalty = 500
    warning('Not starting in 1977 or 1968')
  }
  
  
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
  if(imprinting.switch == TRUE){
    Hm = pars['Hm'] 
  }else{
    Hm = NULL
  }
  if(mass.action == TRUE){
    cc = pars['cc'] 
  }else{
    cc = NULL
  }
  
  
  if(!is.null(Hm)){
    if(Hm > 1){ 
      Hm = -1 
    }} # This will return an infite likelihoos in the line below
  
  if(any(c(alpha.hat, AA.hat, R0.hat, vv.hat, Hm, tau.hat) < 0)) {
    nll.yearly = rep(100000, 10) # Strong penalty for any parameter values that go negative
  }else{
    
    
    if(mass.action == TRUE){
      waifw.in = cc*WAIFW.in + (1-cc)*(WAIFW.in*0+1) # combine some amount of mass action with POLYMOD
    }
    
    # __________________ Use model to generate expected probabilities 
    #___________________ that any case falls in a given birth year ____________________ #  
    sim = simulate.incidence(alpha.hat, AA.hat, R0.hat, vv.hat, tau.hat, Hm, demog = demog, clusters = clusters, start.year = st.yr, end.year = end.yr, WAIFW = WAIFW.in, vaccinate = vaccinate.switch, Fraction.vaccinated = vax.matrix, OAS = OAS.switch, imprinting = imprinting.switch, wm.in = wm.in)
    
    
    
    #___________________ calculate neg log lk ____________________ #  
    fit.years = rownames(dat.in)
    p.hat = sim[fit.years, ]
    nll.yearly = vector('numeric', length(fit.years))
    counter = 0
    for(ii in 1:length(fit.years)){
      if( sum(p.hat[ii, ]) == 0 ){
        #print(ii)
        counter = counter + 1
        nll.yearly[ii] = penalty  # Penalize if no outbreak occurs in the year of interest
      }else{
        nll.yearly[ii] = -dmultinom(x = as.numeric(dat.in[ii,]), prob = as.numeric(p.hat[ii,]), log = TRUE)
      }
      
      
      if(plot == TRUE ){
        cols = rainbow(length(fit.years))
        if(ii == 1){
          plot(0:99, p.hat[ii, ], xlab = 'Age', ylab = 'Fraction of cases', main = paste(round(R0.hat, 2), round(AA.hat, 2), round(alpha.hat, 2), round(Hm, 2)), ylim = c(0, .05), col = cols[ii], type = 'l')
          points(0:99, dat.in[ii, ]/sum(dat.in[ii, ]), col = cols[ii])
        }else{
          lines(0:99, p.hat[ii, ], col = cols[ii])
          points(0:99, dat.in[ii, ]/sum(dat.in[ii, ]), col = cols[ii])
        }
        if(ii == length(fit.years)){
          text(60, .04, paste(sum(nll.yearly)))
          text(60, .035, paste(penalty))
        }
      }
      
      
    }
  }
  
  #list(value = sum(nll.yearly), penalty.count = counter)
  sum(nll.yearly)
}

