# ## This script works up to 2017!
# ## This script assumes vaccination delays imprinting!
# 
# 

#Load year-of-infection intensities
intensities = read.csv('~/Dropbox/R/2017_seasonal_flu/Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2017
#Weight the annual probability of infection by intensity
load('~/Dropbox/R/2017_seasonal_flu/pest.RData')
weighted.attack.rate = p.est*(intensities$Intensity); names(weighted.attack.rate) = 1911:2017
weighted.attack.rate = c(rep(p.est, 1910-1830+1), weighted.attack.rate); names(weighted.attack.rate) = 1830:2017
# 



get.type.weights.AB.vaccination_2 = function(years.out, Countries.out, vax.matrix, type = 1, earliest.birth.year = 1918){
  #Vaccination vector must be a named, country-specific vector for years 1918:2015
  #If multiple countries, inpt a matrix of vaccination vectors, one country per row.

  
  #Get reference tables for each country
  source('~/Dropbox/R/2017_Branching_HA_Imprinting/Clean_CocirculationImport_2017.R')
  #source('Clean_CocirculationImport.R')
  
  birth.years = earliest.birth.year:2017
  infection.years = birth.years
  
  
  
  
  #Initialize master weight matrix
  #Rows - country years
  #Cols - birth years
  weights.master.1 = matrix(NA, nrow = length(Countries.out)*length(years.out), ncol = length(birth.years), dimnames = list(paste(rep(years.out, length(Countries.out)), rep(Countries.out,    each = length(years.out)), sep = ''), rev(birth.years)))
  weights.master.2 = weights.master.naiive = weights.master.vaccine = weights.master.3 = weights.master.1
  
  #Repeat for each country of interest
  for(cc in 1:length(Countries.out)){ 
    
    if(Countries.out[cc] == 'China'){cocirculation.dat = get.cocirculation.ref('China')}
    if(Countries.out[cc] == 'Egypt'){cocirculation.dat = get.cocirculation.ref('Egypt')}
    if(Countries.out[cc] == 'Indonesia'){cocirculation.dat = get.cocirculation.ref('Indonesia')}
    if(Countries.out[cc] == 'Cambodia'){cocirculation.dat = get.cocirculation.ref('Cambodia')}
    if(Countries.out[cc] == 'Vietnam'){cocirculation.dat = get.cocirculation.ref('Vietnam')}
    if(Countries.out[cc] == 'Thailand'){cocirculation.dat = get.cocirculation.ref('Thailand')}
    if(Countries.out[cc] == 'USA'){cocirculation.dat = get.cocirculation.ref('USA')}
    
    early.cocirculation = rbind( c(rep(1, 18), rep(0, 11), rep(1, 59)), #H1
                                 rep(0, 19+11+58), #H2
                                 c(rep(0, 18), rep(1, 11), rep(0, 59)),
                                 c(rep(1, 18), rep(0, 11), rep(1, 59)),
                                 c(rep(0, 18), rep(1, 11), rep(0, 59))) #H3
    
    # add pre-1918 circulation data if necessary
    colnames(early.cocirculation) = 1917:1830
    if(earliest.birth.year < 1918){
      cocirculation.dat = cbind(cocirculation.dat, early.cocirculation[,as.character(1917:earliest.birth.year)])
    }
    
    
    #Extract and re-order relevant data
    #These describe the fraction of circulating viruses in each year that belong to group1 or group2
    type1.dat = (cocirculation.dat['H1', as.character(birth.years)]) 
    type2.dat = (cocirculation.dat['H2', as.character(birth.years)]) 
    type3.dat = (cocirculation.dat['H3', as.character(birth.years)]) 
    
    ## THIS FUNCTION RETURNS A LIST OF:
    # 1. The annual probabilities of infection for a given birth year, from the perspective of a given incidence year
    # 2. The naiive fraction if the birth year is within 12 years of the incidence year
    
    ## Calculate e_ij, prob of any infection in a given year given birth in year i
    #Output a vector of probabilities of being first infected in any year, j given birth in year i
    source('~/Dropbox/R/2017_seasonal_flu/gete_ij_vaccination.R')
    
    ### This function fills in a matrix of exposure weights (equations 2 and 3 in manuscript), for each year of reference
    # This function calls function get.e_ij
    # OUTPUT - A matrix of all probabilites, e_ij with birth years on rows, and years of first infection on columns
    # INPUT - The incidence year gives the year of reference (e.g. H5N1 infections occurring in 2003 would still have naiive people back to 1991, wheras H5N1 infections occurring in 2015 would have naiive individuals going back only to 2003)
    get.infection.weights = function(years.out, type){
      
      vax.vector.full = vax.matrix[,Countries.out[cc]] #This stores the vector contained in the country-specific column
      
      #skips = 2015-incidence.year+1
      H1.mat = matrix(0, nrow = length(years.out), ncol = length(birth.years), dimnames = list((years.out), (birth.years)))
      naiive.mat = H2.mat = H3.mat = H1.mat
      for(jj in 1:length(years.out)){
        for(ii in 1:(years.out[jj]-earliest.birth.year+1)){
          n.inf.years = min(12, years.out[jj]-birth.years[ii])
          n.vax.years = min(10, years.out[jj]-birth.years[ii]+1)
          use.row = which(vax.matrix$year == birth.years[ii])

        #options(warn = 2)
          
          vax.vector.birth_year_specific = rep(0, 13) #Initialize with all 0s
          vax.vector.birth_year_specific[1:(n.vax.years)] = vax.vector.full[use.row:(use.row+n.vax.years-1)] #Fill in first 6 entries with relevant coverage

          inf.probs = get.e_ij.vax(birth.years[ii], years.out[jj], vax.vector = vax.vector.birth_year_specific)
          #vaccine.probs = get.e_ij.vax_causing_imprint(birth.years[ii], years.out[jj], vax.vector = vax.vector.birth_year_specific)$lambda_ij
          
          
          
          
          if(length(inf.probs) == 13){
            #n.factor = sum(inf.probs, vaccine.probs)
            inf.probs = inf.probs/sum(inf.probs)
            #vaccine.probs = vaccine.probs/n.factor
          }
          
          
          
          
          
          H1.mat[jj, ii] = sum(inf.probs*type1.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
          H2.mat[jj, ii] = sum(inf.probs*type2.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
          H3.mat[jj, ii] = sum(inf.probs*type3.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
          #vaccinated.mat[jj, ii] = sum(vaccine.probs)
          naiive.mat[jj, ii] = round(1-sum(inf.probs), digits = 8) #Rounds to the nearest 8 to avoid machine 0 errors
        }
      }
      #return the output in order of 2015:2918
      H1.mat = H1.mat[,as.character(max(birth.years):min(birth.years))]
      H2.mat = H2.mat[,as.character(max(birth.years):min(birth.years))]
      H3.mat = H3.mat[,as.character(max(birth.years):min(birth.years))]
      #vaccinated.mat = vaccinated.mat[,as.character(max(birth.years):min(birth.years))]
      naiive.mat = naiive.mat[,as.character(max(birth.years):min(birth.years))]
      
      
#       if(type == 1){return(H1.mat)}
#       if(type == 2){return(H2.mat)}
#       if(type == 3){return(H3.mat)}
#       if(type == 4){return(naiive.mat)}
return(list(H1.mat = H1.mat, H2.mat = H2.mat, H3.mat = H3.mat, naiive.mat = naiive.mat))
    }
    
    #Fill in the appropriate country-specific rows of the master matrix
wts = get.infection.weights(years.out = years.out)
    weights.master.1[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H1.mat
    weights.master.2[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H2.mat
    weights.master.3[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H3.mat
    weights.master.naiive[((cc-1)*length(years.out))+1:length(years.out), ] = wts$naiive.mat
    #weights.master.vaccine[((cc-1)*length(years.out))+1:length(years.out), ] = wts$vaccinated.mat
  }
  
  if(type == 1){ return(weights.master.1)}
  if(type == 2){ return(weights.master.2)}
  if(type == 3){return(weights.master.3)}
  if(type == 4){ return(weights.master.naiive)}
  if(type == 5){return(list(weights.master.1 = weights.master.1, weights.master.2 = weights.master.2, weights.master.3 = weights.master.3, weights.master.naiive = weights.master.naiive))}
  
}


