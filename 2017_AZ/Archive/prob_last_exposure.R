

q_i = c(.3, .1, .2, .6, .09, .07) # Gives probs of any past exposure to cluster i  (i %in% 1, 2, 3,... 6)

##########################################################
## For loop method
##########################################################
p_x = vector('double', length = length(q_i)+1) # Vector of probabilities that the most recent exposure was to to cluster i. First entry represents probability of no previous exposure, second entry represents probability of most recent exposure to cluster 1, etc.

# General strategy:
#   p_x defines the probability of most recent exposure to cluster k
#   k can take values: {NULL, 1, 2, 3, ... N}, where N is the total number of clusters that have circulated
#   P(K = k) = q_k * prod_{j > k} (1-q_j)          ie: prob exposure to k * prob no exposure to any cluster that circulated after k

p_x[1] = prod(1-q_i) # First entry is the null set (no exposure to any cluster)
p_x[length(p_x)] = q_i[length(q_i)] # Last entry is just the probability of exposure to the last cluster that has circulated
# Middle entries follow the formula above
n.clust = length(q_i)
for(ii in 2:n.clust){
  p_x[ii] = q_i[ii-1]*prod(1-q_i[ii:n.clust])
}

# Check that probabilities sum to 1
if(round(sum(p_x), 15) != 1) warning('Probabilities of last exposure do not sum to 1')



##########################################################
## Vectorized method
##########################################################
# Create binary matrices. Multiply down columns to vectorize the above calculation
last.exposure.indicator = t(rbind(rep(0, length(q_i)), diag(length(q_i)))) # Each column represents a different possible last exposure, starting with null
lack.of.exposure.indicator = t(rbind(rep(1, length(q_i)), (upper.tri(matrix(NA, length(q_i), length(q_i)))))) # Each column represents the clusters an individual must not have been exposed to
lack.of.exposure.indicator[which(last.exposure.indicator + lack.of.exposure.indicator == 0)] = NA

# Fully vectorized using a log transformation
p_x_2 = exp(colSums(log(q_i*last.exposure.indicator + (1-q_i)*lack.of.exposure.indicator), na.rm = TRUE) )

# Use apply to multiply downt the columns
p_x_3 = apply(q_i*last.exposure.indicator + (1-q_i)*lack.of.exposure.indicator, 2, function(xx) prod(xx, na.rm = TRUE))

# Check that probabilities sum to 1
if(round(sum(p_x_2), 15) != 1) warning('Probabilities of last exposure do not sum to 1')
if(round(sum(p_x_3), 15) != 1) warning('Probabilities of last exposure do not sum to 1')

# Check that all calculations give the same result
p_x - p_x_2
p_x - p_x_3
p_x_2 - p_x_3
# Some rounding errors in the faster methods. Let's see what the difference in speed is:


system.time( ## For loop method
  for(ii in 1:10000){
    p_x = vector('double', length = length(q_i)+1) # Vector of probabilities that the most recent exposure was to to cluster i. First entry represents probability of no previous exposure, second entry represents probability of most recent exposure to cluster 1, etc.
    
    # General strategy:
    #   p_x defines the probability of most recent exposure to cluster k
    #   k can take values: {NULL, 1, 2, 3, ... N}, where N is the total number of clusters that have circulated
    #   P(K = k) = q_k * prod_{j > k} (1-q_j)          ie: prob exposure to k * prob no exposure to any cluster that circulated after k
    
    p_x[1] = prod(1-q_i) # First entry is the null set (no exposure to any cluster)
    p_x[length(p_x)] = q_i[length(q_i)] # Last entry is just the probability of exposure to the last cluster that has circulated
    # Middle entries follow the formula above
    n.clust = length(q_i)
    for(ii in 2:n.clust){
      p_x[ii] = q_i[ii-1]*prod(1-q_i[ii:n.clust])
    }
    
    # Check that probabilities sum to 1
    if(round(sum(p_x), 15) != 1) warning('Probabilities of last exposure do not sum to 1')
  }
)



system.time( ## Vectorized method
  for(ii in 1:10000){
    # Fully vectorized using a log transformation
    p_x_2 = exp(colSums(log(q_i*last.exposure.indicator + (1-q_i)*lack.of.exposure.indicator), na.rm = TRUE) )
    
    if(round(sum(p_x_2), 15) != 1) warning('Probabilities of last exposure do not sum to 1')
  }
)



system.time( ## Vectorized method
  for(ii in 1:10000){
    # Use apply to multiply downt the columns
    p_x_3 = apply(q_i*last.exposure.indicator + (1-q_i)*lack.of.exposure.indicator, 2, function(xx) prod(xx, na.rm = TRUE))
    
    # Check that probabilities sum to 1
    if(round(sum(p_x_3), 15) != 1) warning('Probabilities of last exposure do not sum to 1')
  }
)



### For loop is actually fastest, and doesn't give rounding errors so let's use that.

