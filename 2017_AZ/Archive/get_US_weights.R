#Load year-of-infection intensities
intensities = read.csv('Intensitymatser.csv', col.names = c('Year', 'Intensity')); rownames(intensities) = 1911:2015
#Weight the annual probability of infection by intensity
load('pest.RData')
weighted.attack.rate = p.est*(intensities$Intensity); names(weighted.attack.rate) = 1911:2015




get.type.weights.AB.vaccination = function(years.out, Countries.out, vaccination.vector, type = 1){
  #Vaccination vector must be a named, country-specific vector for years 1918:2015
  #If multiple countries, inpt a matrix of vaccination vectors, one country per row.
  
  #Get reference tables for each country
  source('Clean_CocirculationImport.R')
  #source('Clean_CocirculationImport.R')
  
  birth.years = 1918:2015
  infection.years = birth.years
  
  
  
  
  #Initialize master weight matrix
  #Rows - country years
  #Cols - birth years
  weights.master.1 = matrix(NA, nrow = length(Countries.out)*length(years.out), ncol = length(birth.years), dimnames = list(paste(rep(years.out, length(Countries.out)), rep(Countries.out,    each = length(years.out)), sep = ''), rev(birth.years)))
  weights.master.2 = weights.master.naiive = weights.master.3 = weights.master.1
  
  #Repeat for each country of interest
  for(cc in 1:length(Countries.out)){ 
    
    if(Countries.out[cc] == 'China'){cocirculation.dat = get.cocirculation.ref('China')}
    if(Countries.out[cc] == 'Egypt'){cocirculation.dat = get.cocirculation.ref('Egypt')}
    if(Countries.out[cc] == 'Indonesia'){cocirculation.dat = get.cocirculation.ref('Indonesia')}
    if(Countries.out[cc] == 'Cambodia'){cocirculation.dat = get.cocirculation.ref('Cambodia')}
    if(Countries.out[cc] == 'Vietnam'){cocirculation.dat = get.cocirculation.ref('Vietnam')}
    if(Countries.out[cc] == 'Thailand'){cocirculation.dat = get.cocirculation.ref('Thailand')}
    if(Countries.out[cc] == 'USA'){cocirculation.dat = get.cocirculation.ref('USA')}
    
    
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
    source('gete_ij_vaccination.R')
    
    ### This function fills in a matrix of exposure weights (equations 2 and 3 in manuscript), for each year of reference
    # This function calls function get.e_ij
    # OUTPUT - A matrix of all probabilites, e_ij with birth years on rows, and years of first infection on columns
    # INPUT - The incidence year gives the year of reference (e.g. H5N1 infections occurring in 2003 would still have naiive people back to 1991, wheras H5N1 infections occurring in 2015 would have naiive individuals going back only to 2003)
    get.infection.weights = function(years.out, type){
      #skips = 2015-incidence.year+1
      H1.mat = matrix(0, nrow = length(years.out), ncol = length(birth.years), dimnames = list((years.out), (birth.years)))
      naiive.mat = H2.mat = H3.mat = H1.mat
      for(jj in 1:length(years.out)){
        for(ii in 1:(years.out[jj]-1918+1)){
          n.inf.years = min(12, years.out[jj]-birth.years[ii])
          inf.probs = get.e_ij.vax(birth.years[ii], years.out[jj], vax.vector = vaccination.vector[as.character(birth.years[ii]:years.out[jj])])
          if(length(inf.probs) == 13) inf.probs = inf.probs/sum(inf.probs)
          H1.mat[jj, ii] = sum(inf.probs*type1.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
          H2.mat[jj, ii] = sum(inf.probs*type2.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
          H3.mat[jj, ii] = sum(inf.probs*type3.dat[as.character(birth.years[ii:(ii+n.inf.years)])])
          naiive.mat[jj, ii] = round(1-sum(inf.probs), digits = 8) #Rounds to the nearest 8 to avoid machine 0 errors
        }
      }
      #return the output in order of 2015:2918
      H1.mat = H1.mat[,as.character(max(birth.years):min(birth.years))]
      H2.mat = H2.mat[,as.character(max(birth.years):min(birth.years))]
      H3.mat = H3.mat[,as.character(max(birth.years):min(birth.years))]
      naiive.mat = naiive.mat[,as.character(max(birth.years):min(birth.years))]
      
      
#       if(type == 1){return(H1.mat)}
#       if(type == 2){return(H2.mat)}
#       if(type == 3){return(H3.mat)}
#       if(type == 4){return(naiive.mat)}
return(list(H1.mat = H1.mat, H2.mat = H2.mat, H3.mat = H3.mat, naiive.mat = naiive.mat))
    }
    
#     #Fill in the appropriate country-specific rows of the master matrix
#     weights.master.1[((cc-1)*length(years.out))+1:length(years.out), ] = get.infection.weights(years.out = years.out, type = 1)
#     weights.master.2[((cc-1)*length(years.out))+1:length(years.out), ] = get.infection.weights(years.out = years.out, type = 2)
#     weights.master.3[((cc-1)*length(years.out))+1:length(years.out), ] = get.infection.weights(years.out = years.out, type = 3)
#     weights.master.naiive[((cc-1)*length(years.out))+1:length(years.out), ] = get.infection.weights(years.out = years.out, type = 4)

#Fill in the appropriate country-specific rows of the master matrix
wts = get.infection.weights(years.out = years.out) #Get all weight types for the year of interest

weights.master.1[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H1.mat #Fill in the appropriate country-year row of the weights.master.## matrix
weights.master.2[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H2.mat
weights.master.3[((cc-1)*length(years.out))+1:length(years.out), ] = wts$H3.mat
weights.master.naiive[((cc-1)*length(years.out))+1:length(years.out), ] = wts$naiive.mat
  }
  
  if(type == 1){ return(weights.master.1)}
  if(type == 2){ return(weights.master.2)}
  if(type == 3){return(weights.master.3)}
  if(type == 4){ return(weights.master.naiive)}
  if(type == 5){return(list(weights.master.1, weights.master.2, weights.master.3, weights.master.naiive))}
  
}


