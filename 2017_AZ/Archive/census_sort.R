## Census sort funciton
#  Input census data by year and sort into specified age classes
# 
#  INPUTS:
#     counts - vector containing # of people of a given age, from 1 to the max age
#     bin.edges - vector of age class boundaries

sort.birth.years = function(rawdata){
  newdat = as.data.frame(matrix(NA, ncol = 4, nrow = 101, dimnames = list(NULL, c('Year', 'Age', 'BirthYear', 'Population'))))
  newdat$Age = 0:100
  newdat$Year = rep(rawdata[1,1], 101)
  newdat$Population[1:85] = rawdata[1:85,3]
  newdat$Population[86:100] = rep(rawdata[86:88,3], each = 5)/5
  newdat$Population[101] = rawdata[89,3]
  newdat$BirthYear = newdat$Year - newdat$Age
  newdat
}

census.sort = function(counts, bin.edges){
  n.bins = length(bin.edges)-1 #Set # age classes
  age.str = vector('numeric', n.bins) #Create an empty vector with the right number of age classes
  #Sum the number of people in each age class
  for(ii in 1:n.bins){
    valid = seq(bin.edges[ii], bin.edges[ii+1]-1)+1 #Add 1 for 0 term
    age.str[ii] = sum(counts[valid])
  }
  age.str
}


freq.table = function(values, bin.edges){
  n.bins = length(bin.edges)-1
  freq = vector('numeric', n.bins)
  for(ii in 1:n.bins){
    freq[ii] = sum(values >= bin.edges[ii] & values < bin.edges[ii+1])
  }
  freq
}

# Calculate the expected incidence
exp.counts = function(country.incidence, country.census){
  total.counts.by.year = rowSums(country.incidence)
  dem.percents = country.census/rowSums(country.census)
  exp.counts = dem.percents*total.counts.by.year
  list(expected = exp.counts, diff = country.incidence - exp.counts)
}

single.year.sort = function(csv.dat){
  
  ### INDIVIDUAL YEARS
  age.str.years = csv.dat$Both.Sexes.Population
  old = rep(age.str.years[(length(age.str.years)-3):length(age.str.years)], each = 5)/5 #Repeat older years
  age.str.years = age.str.years[-((length(age.str.years)-3):length(age.str.years))] #Remove grouped ends
  age.str.years = c(age.str.years, old); names(age.str.years) = seq(0, 104) #Replace with old
  age.str.years
}
  
####################################
#Split yearly vector up into groups
##-------------------------------------
#Ages.by.year is a vector giving the number of individuals of a given age, 0-104
# -- MUST INCLUDE NAMES SPECIFYING YEAR
#gp.edges specifies the age group divisions (minimum and all upper bounds)
census.groups = function(ages.by.year, gp.edges){
  
  if(is.null(dim(ages.by.year))){   # ------------ #If input is a vector
  census.groups = numeric(length(gp.edges)-1)
  nms = character(length(gp.edges)-1)
  for(ii in 1:(length(gp.edges)-1)){
    valid = gp.edges[ii]:(gp.edges[ii+1]-1); valid
    census.groups[ii] = sum(ages.by.year[valid+1]) #use valid+1 to account for zero
    nms[ii] = paste(gp.edges[ii], '-', gp.edges[ii+1]-1, sep = '')
    
  }    
  }else{                # ---------- # If input is a matrix
    census.groups = matrix(NA, dim(ages.by.year)[1], length(gp.edges)-1)
    nms = character(length(gp.edges)-1)
    for(ii in 1:(length(gp.edges)-1)){
      valid = gp.edges[ii]:(gp.edges[ii+1]-1); valid
      census.groups[,ii] = rowSums(ages.by.year[,valid+1]) #use valid+1 to account for zero
      nms[ii] = paste(gp.edges[ii], '-', gp.edges[ii+1]-1, sep = '')

}
}
colnames(census.groups) = nms
census.groups
}

## --------------------------------------
