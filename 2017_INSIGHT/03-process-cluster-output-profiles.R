#### Process likelihood profile outputs form UCLA's Hoffman2 cluster
rm(list = ls())
setwd('~/Dropbox/R/2018_seasonal_flu/2017_INSIGHT/')
load('processed-data/fitted_models.RData') # Load fitted models
library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)
# Import a long-form text file of all profiles
long_pro = as.data.frame(read.table('cluster_outputs/INSIGHT_profiles.txt'))
#long_pro = as.data.frame(read.table('cluster_outputs/test.txt'))
# Add column names
names(long_pro) = c('grid.val', 'prof.nll', 'modname', 'parname', 'date')


#### OUTPUTS
outfile1 = 'processed-data/INSIGHT_par_ests_CIs.csv'
outfile2 = 'processed-data/INSIGHT_formatted_fits_table.csv'




################################
## Use ggplot to plot profiles, as well as best estimates, and check that the profiles intersect MLE estimates
################################
# Create table of par estimates from each model
parmat = matrix(NA, nrow = length(fits), ncol = length(fits$lk.ATUVG$par)+2, dimnames = list(names(fits), c('nll', 'delAIC', names(fits$lk.ATUVG$par))))
parmat[,'nll'] = (sapply(fits, FUN = function(xx) xx$value))
parmat[names(del.AIC),'delAIC'] = del.AIC
for(pp in names(fits$lk.ATUVG$par)){ ## Fill in par estimates
parmat[,pp] = (sapply(fits, FUN = function(xx) xx$par[pp]))
}


parmat = as.data.frame(parmat)
parmat$modname = gsub(pattern = 'lk.(\\w+)', replacement = "\\1", rownames(parmat))

## Melt the data frame
long_parmat = melt(as.data.frame(parmat), id.vars = c('nll', 'delAIC', 'modname'))


# # Plot profiles
# for(mm in parmat$modname){
# pp = ggplot()+
#   geom_point(data = subset(long_pro, modname == mm), aes(x = grid.val, y = prof.nll))+
#   geom_point(data = subset(long_parmat, modname == mm), aes(x = value, y = nll), color = 'red')+
#   geom_vline(data = subset(long_parmat, modname == mm), aes(xintercept = value), color = 'red')+
#   facet_wrap(~parname)
# ggsave(filename = paste('scratch-figures/checkINSIGHTprofs_', mm, '.pdf', sep = ""), plot = pp)
# }
# #These all look good except that a few are missing





################################
## Calculate 95% CIs for all values in profile
################################
## First, make a list of each model and paramter combination for which there should be a profile
modparlist = long_parmat[!is.na(long_parmat$value),]

####################################
## Calculate 95% profile CIs for each model and parameter
####################################
# Define a function to calculate the likelihood ratio threshold, as a function of the best nll value, and the number of constrained pars (df = 1 if only one par is fixed in the profile)
LR.Threshold = function(NLL_best, df){
  #-2log(LR) is distributed chi2 with df given by the number of fixed parameters in the profile
  #algebraically, the threshold for being in the CI is given by:
  threshold = qchisq(.95, df)/2+NLL_best
  threshold
}

# Define a function to extract the CI bounds from the grid, and the profile values
LR.CI = function(threshold, nll.vec, pars.vec){
  if(any(which(nll.vec > threshold) < which.min(nll.vec))  &  any(which(nll.vec > threshold) > which.min(nll.vec)) ){ #If the minimum is not an end point
    #Find string before and after the min value
    lower = nll.vec[1:which.min(nll.vec)]
    upper = nll.vec[(which.min(nll.vec)-1):length(nll.vec)]
    #Extract low CI from first string, and upper CI from upper string
    CI = c(pars.vec[which.min(abs(lower-threshold))], pars.vec[length(lower)+which.min(abs(upper-threshold))])
  }else{
    #If the first value is the minimum
    if(any(which(nll.vec > threshold) < which.min(nll.vec)) == FALSE){ CI = c(pars.vec[1], pars.vec[which.min(abs(nll.vec-threshold))]) }
    #If the last value is the maximum
    if(any(which(nll.vec > threshold) > which.min(nll.vec)) == FALSE){ CI = c(pars.vec[which.min(abs(nll.vec-threshold))], pars.vec[length(nll.vec)])}
  }
  CI
}


## Extract relevant profile data
LR.CI.wrapper = function(mod.in, par.in){
valid = long_pro[long_pro$modname == mod.in & long_pro$parname == par.in, ]
LR.CI(threshold = LR.Threshold(NLL_best = parmat[paste('lk.', mod.in, sep = ''), 'nll'], df = 1), 
       nll.vec = valid$prof.nll, pars.vec = valid$grid.val)
}

## Apply the wrapper across each row in the modparlist
CIs=mapply(FUN = LR.CI.wrapper, mod.in = modparlist[,'modname'], par.in = modparlist[,'variable'])
rownames(CIs) = c('low', 'high')

## Add to modparlist
modparlist = cbind(modparlist, t(CIs))

#If any parameter-model combinations are missing, the CIs will be output as a list.
#Use the code below to diagnose which combinations are missing, and then fix as needed
#modparlist[which(sapply(CIs, FUN = length) !=2), ]
# Check that all are present
#table(long_pro$modname, long_pro$parname)

# ## If any grids are not wide enough to span the CI, LR.CI.wrapper will output NA
# ## Use this code to diagnose, and expand the tested grid
# modparlist[which(is.na(colSums(CIs))),]


## Once import is good, save outputs
write.csv(modparlist, file = outfile1, row.names = FALSE)


## Reformat as a latex table
pasted = paste(sprintf('%.2f', modparlist$value), ' (', sprintf('%.2f', modparlist$low), '-', sprintf('%.2f', modparlist$high), ')', sep = "")

modparlist$pasted = pasted
# Extract variables of interest
longouts = modparlist[,c('modname', 'variable', 'pasted')]
# Add del.AIC as a variable
das = data.frame(modname = modparlist$modname, variable = rep('delAIC', nrow(modparlist)), pasted = modparlist$delAIC)
#prune duplicates
das = unique(das)
das$pasted = sprintf('%.02f', das$pasted)
longouts = rbind(longouts, das)

# wide format into a table
outs = spread(data = longouts, variable, pasted)
# Sort by del.AIC
outs = outs[order(as.numeric(outs$delAIC)), ]
outs = t(outs)
# put modnames and delAIC on top
outs = outs[c(1,17, 2:16), ]

## Write formatted table
write.csv(outs, outfile2, row.names = TRUE)


## Calculate AIC weights using del.AIC
raw.weights = exp(-.5*del.AIC)
AIC.weights = raw.weights/sum(raw.weights)
AIC.weights

imp.type = sub(pattern = "lk.\\w+?([NGS?])", replacement = "\\1", names(AIC.weights))
imp.type[grep('lk.', imp.type)] = 'none'

## Figure out how much weight each type gets
imp.type.weights = sapply(c('N', 'S', 'G', 'none'), FUN = function(tt){
  valid = AIC.weights[imp.type == tt]
  sum(valid)
})

pie(sort(imp.type.weights))

length(parmat[1,])
paste(parmat[1,], 1:16)
