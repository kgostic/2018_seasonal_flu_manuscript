#############################################
# -----
#####     LIKELIHOOD FUNCTION - Broad immune effects
# -----
#############################################
#calculate the likelihood given demography, exposure and immune effects
# estimate g1, the protective effect of g1 exposure to a g1 virus
# OR fix g1 = 1 to return the null model: demography and exposure only

## ___________________________
## 1. A - run any model with all pars fixed at 0
## ---------------------------


## ___________________________
## 2. AH - Hemagglutinin history only
## ---------------------------
lk.A.H = function(pars, wm, wo, age.spline, xx){
  ##            ---- INPUTS ----
  # wm -- vector of fraction of each birth cohort with matched first exposure to H (i.e. to H1N1 or H2N2)
  # wo -- vector of fraction of each birth cohort with unmatched first H exposure (i.e. H3N2)
  # wn -- vector of fractin of each birth cohort naiive to influenza A
  # p0 = null hypothesis for proportions
  # xx -- the number of cases in each age class
  
  ##          ----- OUTPUTS ----
  # Output the -log likelihood for a given year's data
  
  # 1. Assign paramerets to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  #Ho = pars['Ho']  #Relative susceptibility of the naiive
  
  # 2. calculate p_i as a function of the parameters:
  pp = age.spline*(Hm*wm+wo)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}


## ___________________________
## 2. AN - Neuraminidase history only
## ---------------------------
lk.A.N = function(pars, wm, wo, age.spline, xx){
  ##            ---- INPUTS ----
  # wm -- vector of fraction of each birth cohort with matched first exposure to H (i.e. to H1N1 or H2N2)
  # wo -- vector of fraction of each birth cohort with unmatched first H exposure (i.e. H3N2)
  # wn -- vector of fractin of each birth cohort naiive to influenza A
  # p0 = null hypothesis for proportions
  # xx -- the number of cases in each age class
  
  ##          ----- OUTPUTS ----
  # Output the -log likelihood for a given year's data
  
  # 1. Assign paramerets to be fit
  Nm = pars['Nm']    #Relative susceptibility atrributable to H first exposure
  
  # 2. calculate p_i as a function of the parameters:
  pp = age.spline*(Nm*wm+wo)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}


## ___________________________
## 2. ANH - Hemagglutinin history and NA history
## ---------------------------
lk.A.N.H = function(pars, wmm, wmo, wom, woo, age.spline, xx){
  ##            ---- INPUTS ----
  # wm -- vector of fraction of each birth cohort with matched first exposure to H (i.e. to H1N1 or H2N2)
  # wo -- vector of fraction of each birth cohort with unmatched first H exposure (i.e. H3N2)
  # wn -- vector of fractin of each birth cohort naiive to influenza A
  # p0 = null hypothesis for proportions
  # xx -- the number of cases in each age class
  
  ##          ----- OUTPUTS ----
  # Output the -log likelihood for a given year's data
  
  # 1. Assign paramerets to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  Nm = pars['Nm']
  
  # 2. calculate p_i as a function of the parameters:
  pp = age.spline*(Hm*Nm*wmm + Hm*wmo + Nm*wom + woo)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}



## ___________________________
## 2. AH - Hemagglutinin history, decreasing with age
## ---------------------------
lk.A.H.I = function(pars, wm, wo, age.spline, age.decrease, xx){
  ##            ---- INPUTS ----
  # wm -- vector of fraction of each birth cohort with matched first exposure to H (i.e. to H1N1 or H2N2)
  # wo -- vector of fraction of each birth cohort with unmatched first H exposure (i.e. H3N2)
  # wn -- vector of fractin of each birth cohort naiive to influenza A
  # p0 = null hypothesis for proportions
  # xx -- the number of cases in each age class
  
  ##          ----- OUTPUTS ----
  # Output the -log likelihood for a given year's data
  
  # 1. Assign paramerets to be fit
  Hm = pars['Hm']    #Relative susceptibility atrributable to H first exposure
  #Ho = pars['Ho']  #Relative susceptibility of the naiive
  
  # 2. calculate p_i as a function of the parameters:
  pp = age.spline*((Hm*age.decrease)*wm+wo)
  
  #  3. Likelihood is based on the multinomial density
  if(is.null(dim(xx))){ #DO THIS IF DATA FROM ONE YEAR INPUT AS AN ARRAY
    lk = -dmultinom(xx, size = sum(xx), prob = pp, log = TRUE)
  }else{ #ELSE DO THIS IF MULTI-YEAR DATA INPUT IN A MATRIX
    storage = vector('numeric', dim(xx)[1])
    for(jj in 1:dim(xx)[1]){ #Find the log density for each row (dim 1) and take the negative sum
      storage[jj] = -dmultinom(xx[jj,], size = sum(xx[jj,]), prob = pp[jj,], log = TRUE)
    }
    lk = sum(storage) 
  }
  lk # end function 
}






