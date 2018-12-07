


#-------------------------------------------------------------------------------------------------------------------#
#___________ Define a funciton to calculate the fraction of each age group with last exposure to cluster x _________#
#___________________________________________________________________________________________________________________#
get_last_exposure_probs = function(q_i){
  # Input q_i is a vector of length N whose entries represent the probability of any past exposure to cluster 1, 2, ...N  
  # Output p_x is a vector that gives the fraction of each age group whose most recent exposure was to cluster X, or who have no previous exposures
  p_x = vector('double', length = length(q_i)+1) 
  # First entry represents probability of no previous exposure
  # Second entry represents probability of most recent exposure to cluster 1
  # Third entry represents prob most recent exposure to cluster 2
  # etc.
  
  p_x[1] = prod(1-q_i) # First entry is prob of no exposure to any cluster (null set)
  p_x[length(p_x)] = q_i[length(q_i)] # Last entry is just the probability of exposure to the last cluster that has circulated
  
  # General strategy for calculating p_x:
  #   p_x defines the probability of most recent exposure to cluster k
  #   k can take values: {NULL, 1, 2, 3, ... N}, where N is the total number of clusters that have circulated
  #   P(K = k) = q_k * prod_{j > k} (1-q_j)          ie: prob exposure to k * prob no exposure to any cluster that circulated after k
  
  # Middle entries follow the formula above
  n.clust = length(q_i)
  for(ii in 2:n.clust){
    p_x[ii] = q_i[ii-1]*prod(1-q_i[ii:n.clust])
  }
  
  # Check that probabilities sum to 1
  if(round(sum(p_x), 15) != 1) stop('Probabilities of last exposure do not sum to 1')
  p_x
}











#-------------------------------------------------------------------------------------------------------------------#
#___________ Define a funciton to calculate the fraction of each age group with first exposure to cluster x _________#
#___________________________________________________________________________________________________________________#
get_first_exposure_probs = function(q_i){
  # Input q_i is a vector of length N whose entries represent the probability of any past exposure to cluster 1, 2, ...N  
  # Output p_x is a vector that gives the fraction of each age group whose most recent exposure was to cluster X, or who have no previous exposures
  p_x = vector('double', length = length(q_i)) 
  # First entry represents probability of first exposure to strain 1
  # Second entry represents probability of most recent exposure to cluster 1
  # Third entry represents prob most recent exposure to cluster 2
  # etc.
  
  p_x[1] = q_i[1]  # First entry is prob of no exposure to any cluster (null set)
  
  # General strategy for calculating p_x:
  #   p_x defines the probability of first exposure to cluster k
  #   k can take values: {1, 2, 3, ... N}, where N is the total number of clusters that have circulated
  #   P(K = k) = q_k * prod_{j < k} (1-q_j)          ie: prob exposure to k * prob no exposure to any cluster that circulated after k
  
  # Middle entries follow the formula above
  n.clust = length(q_i)
  for(ii in 2:n.clust){
    p_x[ii] = prod(1-q_i[1:ii-1])*q_i[ii]
  }
  
  # Check that probabilities sum to 1
  if(round(sum(p_x)+prod(1-q_i), 15) != 1) stop('Probabilities of first exposure do not sum to 1')
  p_x
}

