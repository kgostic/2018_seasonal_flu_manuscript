get.e_ij.vax = function(birth.year, incidence.year, vax.vector){
  #vax vector = vaccination fraction from birth year to min(birth.year+12, incidence.year)
  #e_ij = naiive = vector('numeric', 13)
  jjs = birth.year:min(birth.year+12, incidence.year) #Get possible years of infection
  nn = length(jjs) # How many possible years of infection?
  ajs = weighted.attack.rate[as.character(jjs)] #Get corresponding annual attack rates 
  #names(e_ij) = jjs; names(naiive) = jjs
  e_ij = pc_ij = ajs[1]*(1-vax.vector[1])
  if(nn > 1){
  for(ii in 2:nn){
  pc_ij[ii] = pc_ij[ii-1] + ajs[ii]*(1-pc_ij[ii-1])*(1-vax.vector[ii])
  e_ij[ii] = ajs[ii]*(1-pc_ij[ii-1])*(1-vax.vector[ii])
  }
  }
  return(e_ij)
}

# 
# get.e_ij.vax_causing_imprint = function(birth.year, incidence.year, vax.vector){
#   #vax vector = vaccination fraction from birth year to min(birth.year+12, incidence.year)
#   #e_ij = naiive = vector('numeric', 13)
#   jjs = birth.year:min(birth.year+12, incidence.year) #Get possible years of infection
#   nn = length(jjs) # How many possible years of infection?
#   ajs = weighted.attack.rate[as.character(jjs)] #Get corresponding annual attack rates 
#   #names(e_ij) = jjs; names(naiive) = jjs
#   e_ij = e_ij.test = pc_ij = ajs[1]*(1-vax.vector[1])
#   if(nn > 1){
#     for(ii in 2:nn){
#       pc_ij[ii] = pc_ij[ii-1] + ajs[ii]*(1-pc_ij[ii-1])*(1-vax.vector[ii])
#       e_ij[ii] = ajs[ii]*(1-pc_ij[ii-1])*(1-vax.vector[ii])
#       #e_ij.test[ii] = prod((1-ajs[1:(ii-1)])+ajs[1:(ii-1)]*(vax.vector[1:(ii-1)]))*ajs[ii]*(1-vax.vector[ii])
#     }
#   }
#   return(e_ij)
# }
  
  