# ###########################################
# ### This script is designed to calculate likelihood profiles for all models
# ###  fitted to INSIGHT data. Run on the Hoffman2 cluster at UCLA, using an array job
# ##  INPUTS (defined in the .sh bash wrapper, parsed by "optparse" below)
# ##     -outpath   - output file path
# ##     -var       - array run number, used as indexing var in R script
# ###########################################
# 
# ###########################################
# ## Parse inputs from bash wrapper
# ###########################################
# library("optparse") # Load library
# library("parallel")
# option_list = list(
#   # Set up ability to pass an output file path from bash to R
#   make_option(c("-outpath", "--outpath"), type="character", default=NULL, 
#               help="the path to your outpath", metavar="character"),
#   # Set up ability to pass the run number from bash to R
#   make_option(c("-var", "--var"), type = "integer", default = NULL, help="input index number")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser); # Now you have a list called opt with entries opt$var (integer, run number), opt$outpath (character, output file path)
# 
# print(opt$outpath) # Print the output filepath from the cluster for checking

opt = list(outpath = 'cluster_outputs/test.txt', var = 22)

aa = Sys.time()

###########################################
## Load R inputs
###########################################
## Calculate likelihood profiles for all paramters in models fitted to INSIGHT data
source('00-Import_FLU002_-for-multinomial.R')
source('0func-multinomial_likelihood.R')
load('processed-data/fitted_models.RData')
#attach(fits)


#############################
## Write a function wrapper to calculate each profile value in a grid of fixed par values
############################
## Write a function wrapper to calculate one profile point value
## Input pvec, a named vector of best paramter estimates from a fitted likelihood. Should include the paramter to be fixed
##       fixpar, name of the parameter to be fixed
##       fixed.par.vale - value of paramter to be fixed
##       lows -  vector of lower par limits, corresponding to each entry in pvec, including the fixed parameter
##       highs - vector of upper par limits
one_prof_point = function(pvec, fixpar, fixed.par.value, lows, highs, H1.protection, H3.protection){
  drop = which(names(pvec) == fixpar) # Figure out index of fixed par
  # Optimize the likelihood with respect to all free paramters, other than the fixed par
  optim(par = pvec[-drop], fn = profile_function, fixed.par.name = fixpar, fixed.par.value = fixed.par.value, 
        wAV = av.master, wDX = dx.master, wVX = vac.master, 
        wPro.H1 = H1.protection, dat.H1 = H1.master, 
        wPro.H3 = H3.protection, dat.H3 = H3.master, 
        a18.24 = a18.24, a25.31 = a25.31, 
        a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, 
        a67.73 = a67.73, a74.80 = a74.80, a81.90 = a81.90, tested.master = tested.master,
        method = 'L-BFGS-B', lower = lows[-drop], upper = highs[-drop])$value
}


### Set up a second wrapper to apply one_prof_point across a whole grid of fixed parameter values
### INPUTS - 
##     pvec - vector of fitted parameter values from fitted model
##     gridpoints - vector of gridpoints at which you want to fix the fixed paramter
##     fixpar - character, name of paramter value to be fixed
##     lows, highs - vector of lower and upper bounds for free par values, corresponds to parameter order in pvec
##     H1.protection, H3.protection - which protection inputs should we use (NA, HA subtype, HA group)
##     modname - name of fitted model, output to .csv file

### OUTPUTS
##    append outputs to a text file that tracks grid values, profile likelihood vals, model name, fixed paramter and time
prof_grid = function(pvec, fixpar, gridpoints, lows, highs, H1.protection, H3.protection, modname){
  vals = sapply(gridpoints, FUN = one_prof_point, pvec = pvec, fixpar = fixpar, lows = lows, highs = highs, H1.protection = H1.protection, H3.protection = H3.protection)
  # Return a data frame specifying teh model name, par name, prof value, and profile grid points
  ll = length(gridpoints)
  outs = cbind(gridpoints, vals, rep(modname, ll), rep(fixpar, ll), rep(as.character(Sys.time()), ll))
  write(x = t(outs), file = opt$outpath, ncolumns = 5, append = TRUE)
}




#############################
## Write a script that selects one fitted model 
## and calculates a grid of profile values for each 
## parameter in the fitted model
##
## Then we will run this as an array job on Hoffman,
## in which each job runs profiles for a different fitted model.
############################
fitted_model = fits[[opt$var]] # Extract one model for which to run profiles
model_name = gsub(pattern = 'lk.(\\w+)', replacement = '\\1', x = names(fits)[opt$var]) # Extract factor abbreviation from model's name

## Set up lower and upper bounds for estimated free par values in each model. Vector changes depending on which free paramters are included in the model fit
# mod.lows = .001 Always 0.001, hard code below
mod.highs = rep(10, 9) # This defines upper parameter limits for age pars 
if(grepl('[NSG]', model_name)){mod.highs = c(1,1,mod.highs)} # If model includes imprinting, add two upper limits of 1 for pro.H1 and pro.H3
if(grepl('V', model_name)){mod.highs = c(1,1,mod.highs)} # If model includes vaccination, add two upper limits of 1 for rVX.H1 and rVX.H3
if(grepl('U', model_name)){mod.highs = c(10, mod.highs)} # If model includes underlying conditions, add one upper limit of 10 for rDX
if(grepl('T', model_name)){mod.highs = c(10, mod.highs)} # If model includes AV treatment, add one upper limit of 10 for rAV
#rbind(fitted_model$par, mod.highs)

## Set up protection inputs based on model name
modH1pro = modH3pro = 1
if(grepl('N', model_name)){modH1pro = proN1.master; modH3pro = proN2.master} # If the model considers NA subtype-level protection, use NA-specific protection inputs
if(grepl('S', model_name)){modH1pro = proH1.master; modH3pro = proH3.master} # If the model considers HA subtype-level protection, use HA-specific protection inputs
if(grepl('G', model_name)){modH1pro = prog1.master; modH3pro = prog2.master} # If the model considers HA group-level protection, use group-specific protection inputs

## Set up a list of desired grid values at which to calculate profile likelihoods for each free parameter
## Choices are centered around MLE for each parameter
grid.list = list(rAV = seq(.4, 2, by = .005),
                 rDX = seq(.6, 1.5, by = .005),
                 rVX.H1 = seq(.4, 1.3, by = .005),
                 rVX.H3 = seq(.4, 1.3, by = .005),
                 rPro.H1 = seq(.3, 1.3, by = .005),
                 rPro.H3 = seq(.3, 1.3, by = .005),
                 r18.24 = seq(.5, 1.5, by = .005),
                 r25.31 = seq(.5, 1.5, by = .005),
                 r39.45 = seq(.5, 1.5, by = .005),
                 r46.52 = seq(.4, 1.4, by = .005),
                 r53.59 = seq(.4, 1.6, by = .005),
                 r60.66 = seq(.4, 1.6, by = .005),
                 r67.73 = seq(.4, 1.6, by = .005),
                 r74.80 = seq(.1, 1.6, by = .005),
                 r81.90 = seq(.4, 1.6, by = .005))




## Each hoffman array job has 12 threads, which we can use to run profiles in parallel
##    For each array job, calculate profiles for each free parameter in the model.
##    Each parallel job will calculate a profile for one fixed parameter
##       For each parameter, calculate profile values on the grid specified above

## Set up parallel environment
cl = makeCluster(detectCores()-1)
clusterExport(cl, ls()) # Export all variables in workspace
parLapply(cl = cl, X = names(fitted_model$par), 
          fun = function(xx){prof_grid(pvec = fitted_model$par, fixpar = xx, gridpoints = grid.list[[xx]], lows = .001, highs = mod.highs, H1.protection = modH1pro, H3.protection = modH3pro, modname = model_name)})
stopCluster(cl)
#detach(fits)

print(paste(model_name, 'elapsed time =', Sys.time()-aa))


aa = Sys.time()
prof_grid(pvec = fitted_model$par, fixpar = 'r18.24', gridpoints = seq(.5, 1.5, by = .1), lows = .001, highs = mod.highs, H1.protection = modH1pro, H3.protection = modH3pro, modname = model_name)
Sys.time() - aa