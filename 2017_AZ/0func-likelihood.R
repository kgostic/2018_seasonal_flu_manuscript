#######################################
## Define likelihood
######################################
#### INPUTS:
####    pars - named vector of pars to be fit. If not all four pars exist in a given model, code automatically assigns the null value of 1. Parameters not named in "pars" will not be fit.
####    wPro.H1 - n-vector or mxn matrix of the fraction of cases of a given age (n_i), observed in a given country and season (m_i) who are protected against H1 based on childhood imprinting
####    wPro.H3 - same as above, but characterizes protection against seasonal H3N2
####    dat.H1 - n-vector or mxn matrix of observed case counts. Each entry represents the number of H1N1 cases of a given birth year (n_i), observed in a given country and season (m_i).
####    dat.H3 - n-vector or mxn matrix of observed case counts. Each entry represents the number of H3N2 cases of a given birth year (n_i), observed in a given country and season (m_i).
####    a18.24- n-vector or mxn matrix of indicator variables. 1 indicates membership in the 18-24 age class, 0 otherwise
####    All other a##.## inputs are indicators taht follow the same format, but for different age bins.

##### OUTPUTS:
#####    negative log likelihood (scalar)


##### NOTES:  
####    All models tested are nested in this likelihood structure. E.g. to exclude the effects of protection from the model, input wPro = 0, and omit "rPro" from pars.

nll = function(pars, wPro.H1, dat.H1, wPro.H3, dat.H3, a0.4, a5.10, a11.17, a18.24, a25.31, a32.38, a39.45, a46.52, a53.59, a60.66, a67.73, a74.80, a81.90plus){
  # 1. Assign parameters to be fit
  rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
  rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
  b = 1 # Fix relative risk in the baseline group (Ages 0-4) at value 1. Then estimate all others as relative risk. Most should be lower, bounded at 0.
  r5.10 = pars['r5.10'] # Relative risk for 5 to 10 year olds (free paramter to estiamte)
  r11.17 = pars['r11.17'] # Relative risk for 11-17 year olds
  r18.24 = pars['r18.24'] # etc.
  r25.31 = pars['r25.31']
  r32.38 = pars['r32.38']
  r39.45 = pars['r39.45']
  r46.52 = pars['r46.52'] 
  r53.59 = pars['r53.59']
  r60.66 = pars['r60.66'] 
  r67.73 = pars['r67.73'] 
  r74.80 = pars['r74.80'] 
  r81.90p = pars['r81.90p'] 
  
  ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
  age.baseline = b*(a0.4 +
                      r5.10*a5.10 +
                      r11.17*a11.17+
                      r18.24*a18.24+
                      r25.31*a25.31+
                      r32.38*a32.38+ 
                      r39.45*a39.45+ 
                      r46.52*a46.52+ 
                      r53.59*a53.59+ 
                      r60.66*a60.66+ 
                      r67.73*a67.73+ 
                      r74.80*a74.80+ 
                      r81.90p*a81.90plus)
  age.baseline = age.baseline/rowSums(age.baseline) # Normalize so that the fraction of cases predicted in each age group sums to 1 across all age groups
  
  
  # 2. calculate predicted distribution, pp, as a function of the parameters:
  # This step gives the model prediction for H1N1 cases
  pp.H1 = age.baseline * (wPro.H1*rPro.H1+(1-wPro.H1))
  # This step gives the model prediction for H3N2 caeses
  pp.H3 = age.baseline * (wPro.H3*rPro.H3+(1-wPro.H3))
  
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(dat.H1))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
    lk.H1 = -dmultinom(dat.H1, size = sum(dat.H1), prob = pp.H1, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(dat.H1)[1])
    for(jj in 1:dim(dat.H1)[1]){ #Find the neg log density for each row (dim 1) and take the sum
      storage[jj] = -dmultinom(dat.H1[jj,], size = sum(dat.H1[jj,]), prob = pp.H1[jj,], log = TRUE)
    }
    lk.H1 = sum(storage) 
  }
  if(is.null(dim(dat.H3))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
    lk.H3 = -dmultinom(dat.H3, size = sum(dat.H3), prob = pp.H3, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(dat.H3)[1])
    for(jj in 1:dim(dat.H3)[1]){ #Find the neg log density for each row (dim 1) and take the sum
      storage[jj] = -dmultinom(dat.H3[jj,], size = sum(dat.H3[jj,]), prob = pp.H3[jj,], log = TRUE)
    }
    lk.H3 = sum(storage) 
  }
  # Total negative log likelihood is the sum of nll of H3N2 data, and of H1N1 data
  lk.H1+lk.H3 # end function 
}


## Write a function wrapper, so that you only have to input a vector of initial par values and par value limits involved in model comparison
nll.wrapper = function(pars.in, pro.H1, pro.H3, lower.in, upper.in){
  ## Concatenate vector to initialize pars that may or may not be included, with a vector of all age pars (always included)
  pvec = c(pars.in, r5.10 = 1.1, r11.17 = .9, r18.24 = .9, r25.31 = .9, r32.38 = .9, r39.45 = .9, r46.52 = .9, r53.59 = .9, r60.66 =.9, r67.73= .9, r74.80 = .9, r81.90p = .9)
  
  optim(par = pvec, fn = nll, wPro.H1 = pro.H1, wPro.H3 = pro.H3, a0.4 = a0.4, a5.10 = a5.10, a11.17 = a11.17, a18.24 = a18.24, a25.31 = a25.31, a32.38 = a32.38, a39.45 = a39.45, a46.52 = a46.52, a53.59 = a53.59, a60.66 = a60.66, a67.73 = a67.73, a74.80 = a74.80, a81.90plus = a81.90plus, dat.H1 = H1.master, dat.H3 = H3.master, method = 'L-BFGS-B', lower = c(lower.in, rep(.001, 12)), upper = c(upper.in, rep(5, 12)))
}
  
  
  
  
  
### Write a version of the likelihood function that calculates profiles
###  THIS FUNCTION TAKES THE SAME INPUTS AND OUTPUTS AS NLL, ABOVE, EXCEPT...
###       fixed.par.name specifies the name of the fixed parameter
###       fixed.par.value specifies the value of the fixed paramter
###       ## Do not input prof.par as part of the named "pars" vector
###       ## Code will optimize all paramters named in "pars," while fixing the value of the profile.par
profile_func = function(pars, fixed.par.name, fixed.par.value, wPro.H1, dat.H1, wPro.H3, dat.H3, a0.4, a5.10, a11.17, a18.24, a25.31, a32.38, a39.45, a46.52, a53.59, a60.66, a67.73, a74.80, a81.90plus){
    # 1. Assign parameters to be fit (all age paramters, and those named in pars)
    rPro.H1 = ifelse(is.na(pars['rPro.H1']), 1, pars['rPro.H1'])# Relative risk given imprinting protection
    rPro.H3 = ifelse(is.na(pars['rPro.H3']), 1, pars['rPro.H3'])# Relative risk given imprinting protection
    b = 1 # Fix relative risk in the baseline group (Ages 0-4) at value 1. Then estimate all others as relative risk. Most should be lower, bounded at 0.
    r5.10 = pars['r5.10'] # Relative risk for 5 to 10 year olds (free paramter to estiamte)
    r11.17 = pars['r11.17'] # Relative risk for 11-17 year olds
    r18.24 = pars['r18.24'] # etc.
    r25.31 = pars['r25.31']
    r32.38 = pars['r32.38']
    r39.45 = pars['r39.45']
    r46.52 = pars['r46.52'] 
    r53.59 = pars['r53.59']
    r60.66 = pars['r60.66'] 
    r67.73 = pars['r67.73'] 
    r74.80 = pars['r74.80'] 
    r81.90p = pars['r81.90p'] 
    
    #2. If the parameter of interest is meant to be fixed, replace its value with fixed.par.val
    #   Else, keep the value the same
    rPro.H1 = ifelse(fixed.par.name == 'rPro.H1', fixed.par.value, rPro.H1)
    rPro.H3 = ifelse(fixed.par.name == 'rPro.H3', fixed.par.value, rPro.H3)
    r5.10 = ifelse(fixed.par.name == 'r5.10', fixed.par.value, r5.10)
    r11.17 = ifelse(fixed.par.name == 'r11.17', fixed.par.value, r11.17)
    r18.24 = ifelse(fixed.par.name == 'r18.24', fixed.par.value, r18.24)
    r25.31 = ifelse(fixed.par.name == 'r25.31', fixed.par.value, r25.31)
    r32.38 = ifelse(fixed.par.name == 'r32.38', fixed.par.value, r32.38)
    r39.45 = ifelse(fixed.par.name == 'r39.45', fixed.par.value, r39.45)
    r46.52 = ifelse(fixed.par.name == 'r46.52', fixed.par.value, r46.52)
    r53.59 = ifelse(fixed.par.name == 'r53.59', fixed.par.value, r53.59)
    r60.66 = ifelse(fixed.par.name == 'r60.66', fixed.par.value, r60.66)
    r67.73 = ifelse(fixed.par.name == 'r67.73', fixed.par.value, r67.73)
    r74.80 = ifelse(fixed.par.name == 'r74.80', fixed.par.value, r74.80)
    r81.90p = ifelse(fixed.par.name == 'r81.90p', fixed.par.value, r81.90p)
    
    ## Age-specific baseline prediction takes the same form for H1N1 and H3N2. Attempt to explain residual, subtype-specific differences through differences in imprinting history, etc. below.
    age.baseline = b*(a0.4 +
                        r5.10*a5.10 +
                        r11.17*a11.17+
                        r18.24*a18.24+
                        r25.31*a25.31+
                        r32.38*a32.38+ 
                        r39.45*a39.45+ 
                        r46.52*a46.52+ 
                        r53.59*a53.59+ 
                        r60.66*a60.66+ 
                        r67.73*a67.73+ 
                        r74.80*a74.80+ 
                        r81.90p*a81.90plus)
    age.baseline = age.baseline/rowSums(age.baseline) # Normalize so that the fraction of cases predicted in each age group sums to 1 across all age groups
    
    
    # 2. calculate predicted distribution, pp, as a function of the parameters:
    # This step gives the model prediction for H1N1 cases
    pp.H1 = age.baseline * (wPro.H1*rPro.H1+(1-wPro.H1))
    # This step gives the model prediction for H3N2 caeses
    pp.H3 = age.baseline * (wPro.H3*rPro.H3+(1-wPro.H3))
    
    
    #  3. Likelihood is based on the multinomial density
    if(is.null(dim(dat.H1))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
      lk.H1 = -dmultinom(dat.H1, size = sum(dat.H1), prob = pp.H1, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
    }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
      storage = vector('numeric', dim(dat.H1)[1])
      for(jj in 1:dim(dat.H1)[1]){ #Find the neg log density for each row (dim 1) and take the sum
        storage[jj] = -dmultinom(dat.H1[jj,], size = sum(dat.H1[jj,]), prob = pp.H1[jj,], log = TRUE)
      }
      lk.H1 = sum(storage) 
    }
    if(is.null(dim(dat.H3))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS A VECTOR
      lk.H3 = -dmultinom(dat.H3, size = sum(dat.H3), prob = pp.H3, log = TRUE) #This line returns the log multinomial density of the observed data, with expected probabilities governed by model predictions.
    }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
      storage = vector('numeric', dim(dat.H3)[1])
      for(jj in 1:dim(dat.H3)[1]){ #Find the neg log density for each row (dim 1) and take the sum
        storage[jj] = -dmultinom(dat.H3[jj,], size = sum(dat.H3[jj,]), prob = pp.H3[jj,], log = TRUE)
      }
      lk.H3 = sum(storage) 
    }
    # Total negative log likelihood is the sum of nll of H3N2 data, and of H1N1 data
    lk.H1+lk.H3 # end function 
}
