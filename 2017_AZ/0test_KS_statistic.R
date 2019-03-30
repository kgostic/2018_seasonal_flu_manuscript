## We aim to use the Kolmogorov-Smirnov statistic (KS statistic) to compare age distributions of infection in seasons with different antigenic advance.
##  Sample sizes are not equal between seasons, so to compare KS statistic values, we need to the variance based on sample sizes.
##  In a one-sample KS test, if n is the sample size of the data characterizing the observed distribution, and K is the KS statistic, then sqrt(n)*K converges to the Kolmogorv Distribution. The factor sqrt(n) standardizes the variance.
##  In the two-sample KS test, if n and m and the sample sizes of each data set, then sqrt((m*n)/(m+n))*K should coverge to the Kolmogorov distribution. The factor sqrt((m*n)/(m+n)) standardizes the variance.

### GOALS OF THIS SCRIPT:
## Simulate data from arbitrarily chosen underlying distributions, and with known sample size to answer the following questions:
# 1. Does the statistic output by R's built-in ks.test() represent the raw value, K, or the sandardizes value, sqrt(n)*K?
# 2. Verify that standardized K values follow a consistent distribution, regardless of sample size.



#####################################
##  One-sample test. 
#####################################
## Define a function to simulate data and calculate K three ways:
##    --> Define your own method to calculate K. Output the raw value.
##    --> Output the standardized K value (K*sqrt(n))
##    --> Output the value given by R's built-in function
ksman = function(ss){
  d1 = sort(runif(ss))
  nn = length(d1)
  ## Find the ecdf of d1, d2
  Fx = sapply(d1, FUN = function(xx) sum(d1<=xx)/nn)
  Fx_exp = sapply(d1, pnorm)
  c('Routput' = ks.test(d1, y = 'pexp')$statistic, 'raw' = max(abs(Fx-Fx_exp)), 'standardized' = max(abs(Fx-Fx_exp))*sqrt(nn))
}

### Does R output the standardized value or the raw value?
ksman(100)
ksman(10)
ksman(1000)
#########################################################################################################
############### R OUTPUT IS RAW, NOT STANDARDIZED
#########################################################################################################

## Use the function to simulate data using rnorm, and find the KS statistic comparing simulated data to a true normal distribution
par(mfrow = c(2,1))
plot(density(replicate(n = 1000, ksman(10))[2,]), main = 'Raw K value\nn=10') # Plot distribution of KS outputs, raw value, sample size = 10
xls = par('usr')[1:2]
plot(density(replicate(n = 1000, ksman(100))[2,]), main = 'n=100', xlim = xls) # Sample size = 10
#########################################################################################################
############### THESE DISTRIBUTIONS ARE DIFFERENT BECAUSE THE VARIANCE IS NOT STANDARDIZED.
#########################################################################################################


## Repeat using standardized values
par(mfrow = c(2,1))
plot(density(replicate(n = 1000, ksman(10))[3,]), main = 'Standardized K value, K*sqrt(n)\nn=10') # Plot distribution of KS outputs, raw value, sample size = 10
xls = par('usr')[1:2]
plot(density(replicate(n = 1000, ksman(100))[3,]), main = 'n=100', xlim = xls) # Sample size = 10
#########################################################################################################
############### THESE DISTRIBUTIONS ARE THE SAME BECAUSE THE VARIANCE IS NOW STANDARDIZED
#########################################################################################################











##########################################################################
##  REPEAT CHECK FOR 2-SAMPLE TEST, WHERE EACH SAMPLE HAS THE SAME SAMPLE SIZE
##########################################################################
ksman = function(ss){
  d1 = sort(rnorm(ss)) # draw 2 samples of size ss
  d2 = sort(rnorm(ss))
  nn = length(d1)
  ## Find the ecdf of d1, d2
  Fx1 = sapply(d1, FUN = function(xx) sum(d1<=xx)/nn)
  Fx2 = sapply(d2, FUN = function(xx) sum(d2<=xx)/nn)
  ## Combine into a single vector
  cc = rbind(c(d1, d2),
             c(Fx1, rep(NA,nn)),
             c(rep(NA, nn), Fx2))
  ord = order(cc[1,]);   cc = cc[,ord] # Re-order based on sequential d1 and d2 values
  valid = which(is.na(cc[,1])); cc[valid,1] = 0 # ecdf of the row with na in the first column is 0
  ## Fill in NAs with value from the previous column
  for(ii in 2:(nn*2)){
    if(is.na(cc[2,ii])){ cc[2,ii] = cc[2,ii-1] }
    if(is.na(cc[3,ii])){ cc[3,ii] = cc[3,ii-1] }
  }
  # plot(cc[1,], cc[2,], pch = 16, col = 'gray')
  # points(cc[1,], cc[3,], col = 'red')
  c('R' = ks.test(d1, d2)$statistic, 
    'raw' = max(abs(cc[2,]-cc[3,])),
     'standardized' = max(abs(cc[2,]-cc[3,]))*sqrt((nn^2)/(nn+nn)))
}

ksman(100)


test100 = as.data.frame(t(replicate(1000, ksman(100))))
test500 = as.data.frame(t(replicate(1000, ksman(500))))
test1000 = as.data.frame(t(replicate(1000, ksman(1000))))

pdat = rbind(test100, test500, test1000)
pdat$ss = rep(c(100, 500, 1000), each = 1000)

## Plot raw KS stats
ggplot()+
  geom_density(data = pdat, aes(x = raw, fill = factor(ss)), alpha = .7)

## Plot standardized KS stats
ggplot()+
  geom_density(data = pdat, aes(x = standardized, fill = factor(ss)), alpha = .7)
## multiplying by sqrt((nn^2)/(nn+nn)) standardizes the density
##############################################################################################################
## SAME CONCLUSIONS - R output is not standardized, but multiplying by sqrt((m*n)/(m+n)) standardizes the variance
###############################################################################################################






##########################################################################
##  REPEAT CHECK FOR 2-SAMPLE TEST, UNEQUAL SAMPLE SIZE
##     Note - results are robust to changes in the distribution used to simulate data
##########################################################################
## Two-sample, equal sample sizes -- appears not to standardize variance of stat
ksman = function(ss1, ss2 = 1000){
  d1 = sort(rnorm(ss1)) # draw 2 samples of size ss
  d2 = sort(rnorm(ss2, mean = 3, sd = 2))
  nn = length(d1)
  mm = length(d2)
  ## Find the ecdf of d1, d2
  Fx1 = sapply(d1, FUN = function(xx) sum(d1<=xx)/nn)
  Fx2 = sapply(d2, FUN = function(xx) sum(d2<=xx)/mm)
  ## Combine into a single vector
  cc = rbind(c(d1, d2),
             c(Fx1, rep(NA,mm)),
             c(rep(NA, nn), Fx2))
  ord = order(cc[1,]);   cc = cc[,ord] # Re-order based on sequential d1 and d2 values
  valid = which(is.na(cc[,1])); cc[valid,1] = 0 # ecdf of the row with na in the first column is 0
  ## Fill in NAs with value from the previous column
  for(ii in 2:(nn+mm)){
    if(is.na(cc[2,ii])){ cc[2,ii] = cc[2,ii-1] }
    if(is.na(cc[3,ii])){ cc[3,ii] = cc[3,ii-1] }
  }
  # plot(cc[1,], cc[2,], pch = 16, col = 'gray')
  # points(cc[1,], cc[3,], col = 'red')
  c('R' = ks.test(d1, d2)$statistic, 
    'raw' = max(abs(cc[2,]-cc[3,])),
    'standardized' = max(abs(cc[2,]-cc[3,]))*sqrt((nn*mm)/(nn+mm)))
}

ksman(100)


test100 = as.data.frame(t(replicate(1000, ksman(100, 500))))
test500 = as.data.frame(t(replicate(1000, ksman(500, 700))))
test1000 = as.data.frame(t(replicate(1000, ksman(1000, 300))))

pdat = rbind(test100, test500, test1000)
pdat$ss = rep(c(100, 500, 1000), each = 1000)

## Plot raw KS stats
ggplot()+
  geom_density(data = pdat, aes(x = raw, fill = factor(ss)), alpha = .7)

## Plot standardized KS stats
ggplot()+
  geom_density(data = pdat, aes(x = standardized, fill = factor(ss)), alpha = .7)
## multiplying by sqrt((nn^2)/(nn+nn)) standardizes the density
##############################################################################################################
## SAME CONCLUSIONS - R output is not standardized, but multiplying by sqrt((m*n)/(m+n)) standardizes the variance
###############################################################################################################


## Conclusion:
## Output the KS statistic from R's built-in ks.test() function.
## Multiply by sqrt((m*n)/(m+n)) to standardize the variance, where n is sample size of vector 1, and m is sample size of vector 2
## KS stats will follow the same distribution, regardless of sample size