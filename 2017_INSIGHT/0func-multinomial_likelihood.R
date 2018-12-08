#######################################
## Define likelihood
######################################
#### INPUTS:
####    pars - named vector of pars to be fit. If a parameter is not named in this input vector, code automatically assigns the null value of 1, excluding the factor of interest from the fitted model. Because the models are nested, assigning a parameter value of 1 (e.g. for relative risk given antivial treatment), means that antivial treatment has no effect of overall predicted risk. Thus, the model reduces algebraically to a simpler predictive model that does not consider antiviral treatment effects. Thus, we can use the single likelihood function below to fit all combinations of tested factors. Only factors whose free paramteres are named in the input vector will be included in the overall model fit.
####    tesed.master - matrix of the number of cases tested for influenza in a given country-season, and of a given age. Vectors from each country-season represent the null prediction for that country-season; if all patients were at equal risk from a given subtype, then the observed age distribution of cases of H1N1 or H3N2 would be proportional to the distribution of tested cases.
####    wAV - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_j) who used antivirals
####    wDX - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_j) who had underlying conditions
####    wVX - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_j) who were vaccinated in the previous 12 months
####    wPro - n-vector or mxn matrix of "weights" describing the fraction of cases of a given age (n_i), observed in a given country and season (m_j) who are protected against the challenge strain based on childhood imprinting
####    xx - n-vector or mxn matrix of observed case counts. Each entry represents the number of cases of a given age (n_i), observed in a given country and season (m_j).

##### OUTPUTS:
#####    negative log likelihood (scalar)


##### NOTES:  
####    All models tested are nested in this likelihood structure. E.g. to exclude the effects of protection from the model, input wPro = 0, and omit "rPro" from pars.

nll = function(pars, wAV, wDX, wVX, wPro.H1, dat.H1, wPro.H3, dat.H3, a18.24, a25.31, a32.38, a39.45, a46.52, a53.59, a60.66, a67.73, a74.80, a81.90, tested.master){
  # 1. Assign parameters to be fit
  rAV = ifelse(is.na(pars['rAV']), 1, pars['rAV']) # Relative risk given antiviral use
  rDX = ifelse(is.na(pars['rDX']), 1, pars['rDX']) # Relative risk given underlying symptoms
  rVX.H1 = ifelse(is.na(pars['rVX.H1']), 1, pars['rVX.H1']) # Relative risk given vaccination
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rVX.H3 = ifelse(is.na(pars['rVX.H3']), 1, pars['rVX.H3']) # Relative risk given vaccination
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  r18.24 = pars['r18.24'] # Baseline expectation for 18-24 year olds
  r25.31 = pars['r25.31']
  b =1 # Baseline risk (ages 32-38) are all fixed to 1. All other r##.## paramters represent relative risk, relative to baseline. The choice of baseline risk group is arbitrary. I chose the 32-38 age group, because they are usually at or near the highest risk, and it is somewhat easier to interpret relative risk estimates if they are all above or below 1.
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73'] 
  r74.80 = pars['r74.80'] 
  r81.90 = pars['r81.90'] 
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.risk = b*(r18.24*a18.24+
                  r25.31*a25.31+ 
                  a32.38+ 
                  r39.45*a39.45+ 
                  r46.52*a46.52+ 
                  r53.59*a53.59+ 
                  r60.66*a60.66+ 
                  r67.73*a67.73+ 
                  r74.80*a74.80+ 
                  r81.90*a81.90)
  age.risk = age.risk/rowSums(age.risk) # Normalize so that age specific risk prediction sums to 1 across all age groups.
  
  
  # 2. calculate predicted distribution, pp, as a linear combination of age-specific risk, and other model factors.
  # Risk prediction is a function of the estimated free parameters:
  # This step gives the model predicted age distributions for H1N1 and H3N2
  # No need to normalize the preditions so that probabilities of case observation sum to 1 across age groups. The dmultinom function (see step 3) will normalize predicted probabilities automatically.
  pp.H1 = tested.master * age.risk * (wAV*rAV+(1-wAV))* (wDX*rDX+(1-wDX))* (wVX*rVX.H1+(1-wVX))* (wPro.H1*rPro.H1+(1-wPro.H1))
  
  pp.H3 = tested.master * age.risk * (wAV*rAV+(1-wAV))* (wDX*rDX+(1-wDX))* (wVX*rVX.H3+(1-wVX))* (wPro.H3*rPro.H3+(1-wPro.H3))
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(dat.H1))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
    lk.H1 = -dmultinom(dat.H1, size = sum(dat.H1), prob = pp.H1, log = TRUE) #This line returns the log multinomial density of the observed H1N1 data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(dat.H1)[1])
    for(jj in 1:dim(dat.H1)[1]){ #Find the neg log density for each row (dim 1) and take the sum
      storage[jj] = -dmultinom(dat.H1[jj,], size = sum(dat.H1[jj,]), prob = pp.H1[jj,], log = TRUE)
    }
    lk.H1 = sum(storage) 
  }
  if(is.null(dim(dat.H3))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
    lk.H3 = -dmultinom(dat.H3, size = sum(dat.H3), prob = pp.H3, log = TRUE) #This line returns the log multinomial density of the observed H3N2 data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(dat.H3)[1])
    for(jj in 1:dim(dat.H3)[1]){ #Find the neg log density for each row (dim 1) and take the sum
      storage[jj] = -dmultinom(dat.H3[jj,], size = sum(dat.H3[jj,]), prob = pp.H3[jj,], log = TRUE)
    }
    lk.H3 = sum(storage) 
  }
  ## Sum the negative log likelihoods of the H1N1 and H3N2 data to get the full nll, and return
    lk.H1+lk.H3
# end function 
}

# Test function
nll(pars = c(rAV = .5, rDX = uu, rVX.H1 = .5, rPro.H1 = .5, rVX.H3 = .5, rPro.H3 = .5, r18.24 = .9, r25.31 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73= .9, r74.80 = .9, r81.90 = .9), wAV = av.master, wDX = dx.master, wVX = vac.master, wPro.H1 = proH1.master, wPro.H3 = proH3.master, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90 = a81.90, dat.H1 = H1.master, dat.H3 = H3.master, tested.master = tested.master/rowSums(tested.master))





## Write a function wrapper, so that you can optimize the likelihood, but only have to input a vector of initial par values and par value limits involved in model comparison
nll.wrapper = function(pars.in, pro.H1, pro.H3, lower.in, upper.in){
  ## Concatenate vector to initialize pars that may or may not be included, with a vector of all age pars (always included)
  pvec = c(pars.in, r18.24 = .9, r25.31 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73= .9, r74.80 = .9, r81.90 = .9)
  
  optim(par = pvec, fn = nll, wAV = av.master, wDX = dx.master, wVX = vac.master, wPro.H1 = pro.H1, wPro.H3 = pro.H3, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90 = a81.90, dat.H1 = H1.master, dat.H3 = H3.master, tested.master = tested.master/rowSums(tested.master), method = 'L-BFGS-B', lower = c(lower.in, rep(.001, 9)), upper = c(upper.in, rep(5, 9)))
}


# ## Write another function wrapper so that you can input a vector of parameters and output model predictions for plotting
# pred.wrapper = function(pars.in, pro.H1, pro.H3){
#   ## Concatenate vector to initialize pars that may or may not be included, with a vector of all age pars (always included)
#   pvec = c(rAV = 1, rDX = 1, rVX.H1 = 1, rVX.H3 = 1, rPro.H1 = 1, rPro.H3 = 1, r18.24 = 1, r25.31 = 1, r39.45 = 1, r46.52 = 1, r53.39 = 1, r60.66 = 1, r67.73 = 1, r74.80 = 1, r81.90 = 1) ## defaults
#   ## overwrite input values
#   for(pp in names(pars.in)){
#     pvec[pp] = pars.in[pp]
#   }
#   
#   nll(wAV = av.master, wDX = dx.master, wVX = vac.master, wPro.H1 = pro.H1, wPro.H3 = pro.H3, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90 = a81.90, dat.H1 = H1.master, dat.H3 = H3.master, tested.master = tested.master/rowSums(tested.master), pars = pvec, rtrn = 'prediction')
# }
