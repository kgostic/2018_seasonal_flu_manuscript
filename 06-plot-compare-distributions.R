rm(list = ls())
library(reshape)
setwd('~/Dropbox/R/2018_seasonal_flu/2017_INSIGHT/')
source('00-Import_FLU002_-for-multinomial.R')
setwd('../')

## Plot observed age distributions of infection using the INSIGHT 002 data set
## Repeat using the AZ data set
## MAIN TEXT PLOTS:
##  - Age dists from INSIGHT, overall, and panels showing all seasons and countries with cocirculation
##  - Age dists from AZ, overall and panels with all seasons with cocirculation
## SUPPLEMENTARY TEXT PLOTS
##  - Lower smoothing paramters
##  - Higher smoothing paramters
##  - All countries or all seasons

## Set smoothing paramter for smoothing splines
spar.in = .8
spar.low = .6
spar.high = 1

#### OUTPUTS
outfile1 = 'figures/age-dists-by-country-and-season-INSIGHT.pdf' ## Main text, Fig. 3
outfile2 = 'figures/age-dists-INSIGHT-hiSmoothPar.pdf' ## Supplement
outfile3 = 'figures/age-dists-INSIGHT-lowSmoothPar.pdf' ## Supplement
outfile4 = 'figures/age-dists-INSIGHT-all-countries.pdf' ## Supplement
outfile5 = 'figures/age-dists-INSIGHT-all-seasons.pdf' ## Supplement
outfile6 = 'figures/age-dists-ARIZONA-all-seasons.pdf' ## Supplement
outfile7 = 'figures/age-dists-by-season-ARIZONA-ages0-90.pdf' ## Supplement
outfile8 = 'figures/age-dists-by-season-ARIZONA-ages18-90.pdf' ## Main text, Fig. 2
outfile9 = 'figures/age-dists-by-season-ARIZONA-hilowsmoothpar.pdf' ## Supplement






#################################
## INSIGHT plots
#################################
## Define a function to quickly generate a fraction table for single year of age from 18:90
agerange = 18:90
age.tab = function(agevec){sapply(agerange, FUN = function(aa) sum(agevec == aa))}

## Drop data from the 2009 pandemic
dat.002 = subset(dat.002, !season %in% c('NH.09.10', 'SH.10'))

## Plot only countries and subtypes for which at least 50 cases were observed of each subtype
tb = (table(dat.002$country, dat.002$flutype)) # Get counts for each country, H1N1 = 1, H3N2 = 2
include = rownames(tb)[which(tb[,1]>=50 & tb[,2]>=50)]
## Extract data from countries to be included
cplotdat = subset(dat.002, country %in% include)



### Extract data to plot season-specific distributions
## Plot only seasons for which at least 50 cases were observed of each subtype
tb = (table(dat.002$season, dat.002$flutype)) # Get counts for each season, H1N1 = 1, H3N2 = 2
include = rownames(tb)[which(tb[,1]>=50 & tb[,2]>=50)]
## Extract data from countries to be included
splotdat = subset(dat.002, season %in% include)


## Write a function to output transparent colors
## Credit to: http://www.dataanalytics.org.uk/Data%20Analysis/TipsAndTricks/TTR-20150531.htm
tns <- function(colname, percent = 70) {
  rgb.val <- col2rgb(colname)
  ## Make new color using input color as base and alpha set by transparency
  tnscol <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, alpha = (100-percent)*255/100)
  ## return
  invisible(tnscol)
}



#################################
## Main text INSIGHT plot (Fig. 3)
#################################
{
  pdf(outfile1, height = 9)
  par(mar = c(1,3,2.5,0))
  layout(matrix(c(1,8,2,  8,3,4,  5,6,7), byrow = T, nrow = 3))
  par(mfrow = c(4,3))
  ## Overall
  valid1 = subset(dat.002,flutype == "1" & season!='NH.09.10' & season!='SH.10')
  valid2 = subset(dat.002,flutype == "2" & season!='NH.09.10' & season!='SH.10')
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'), xpd = NA)
  mtext(text = 'age', side = 1, line = 1.9, cex = .8)
  mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
  mtext(text = 'Overall', side = 3, line = .25, font = 2)
  #text(x = 75, y = .08, paste( 'H1N1, n = ', nrow(valid1), '\nH3N2, n = ', nrow(valid2)), cex = .9)
  legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
  ## add smoothed density
  ss = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
  lines(ss$x, ss$y, col = 'dodgerblue')
  points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
  ss = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
  lines(ss$x, ss$y, col = 'firebrick')
  
  
  ### Plot season-specific
  par(mar = c(1.5,3,2,-.5)+1.5)
  ## Country-specific plots
  ssns = unique(splotdat$season)
  for(ss in ssns){
    valid1 = subset(splotdat, season == ss & flutype == "1")
    valid2 = subset(splotdat, season == ss & flutype == "2")
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste(ss), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
      lines(zz$x, zz$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
      zz = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
      lines(zz$x, zz$y, col = 'firebrick')
    }
  }
  
  
  ## Country-specific plots
  ccs = unique(cplotdat$country)
  for(cc in ccs){
    valid1 = subset(cplotdat, country == cc & flutype == "1")
    valid2 = subset(cplotdat, country == cc & flutype == "2")
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text =paste(cc), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      ss = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
      lines(ss$x, ss$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
      ss = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
      lines(ss$x, ss$y, col = 'firebrick')
    }
  }
  dev.off()
}



#################################
## Repeat with high smoothing parameters
#################################
{
  pdf(outfile2, height = 9)
  par(mar = c(1,3,2.5,0))
  layout(matrix(c(1,8,2,  8,3,4,  5,6,7), byrow = T, nrow = 3))
  par(mfrow = c(4,3))
  ## Overall
  valid1 = subset(dat.002,flutype == "1")
  valid2 = subset(dat.002,flutype == "2")
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue', percent = 70), xpd = NA)
  mtext(text = 'age', side = 1, line = 1.9, cex = .8)
  mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
  mtext(text = 'Overall', side = 3, line = .25, font = 2)
  #text(x = 75, y = .08, paste( 'H1N1, n = ', nrow(valid1), '\nH3N2, n = ', nrow(valid2)), cex = .9)
  legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
  ## add smoothed density
  ss1 = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.high) ## High smoothpar
  lines(ss1$x, ss1$y, col = 'dodgerblue')
  points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick', percent = 70))
  ss1 = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.high) ## High smoothpar
  lines(ss1$x, ss1$y, col = 'firebrick')

  
  ### Plot season-specific
  par(mar = c(1.5,3,2,-.5)+1.5)
  ## Country-specific plots
  ssns = unique(splotdat$season)
  for(ss in ssns){
    valid1 = subset(splotdat, season == ss & flutype == "1")
    valid2 = subset(splotdat, season == ss & flutype == "2")
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue', percent = 70))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste(ss), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz1 = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.high)
      lines(zz1$x, zz1$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick', percent = 70))
      zz1 = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.high)
      lines(zz1$x, zz1$y, col = 'firebrick')
    }
  }
  
  
  ## Country-specific plots
  ccs = unique(cplotdat$country)
  for(cc in ccs){
    valid1 = subset(cplotdat, country == cc & flutype == "1")
    valid2 = subset(cplotdat, country == cc & flutype == "2")
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text =paste(cc), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz1 = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.high)
      lines(zz1$x, zz1$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick', percent = 70))
      zz1 = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.high)
      lines(zz1$x, zz1$y, col = 'firebrick')
    }
  }
  dev.off()
}



#################################
## Repeat with low smoothing parameters
#################################
{
  pdf(outfile3, height = 9)
  par(mar = c(1,3,2.5,0))
  layout(matrix(c(1,8,2,  8,3,4,  5,6,7), byrow = T, nrow = 3))
  par(mfrow = c(4,3))
  ## Overall
  valid1 = subset(dat.002,flutype == "1")
  valid2 = subset(dat.002,flutype == "2")
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue', percent = 70), xpd = NA)
  mtext(text = 'age', side = 1, line = 1.9, cex = .8)
  mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
  mtext(text = 'Overall', side = 3, line = .25, font = 2)
  #text(x = 75, y = .08, paste( 'H1N1, n = ', nrow(valid1), '\nH3N2, n = ', nrow(valid2)), cex = .9)
  legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
  ## add smoothed density
  ss1 = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.low) ## Low smoothpar
  lines(ss1$x, ss1$y, col = 'dodgerblue')
  points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick', percent = 70))
  ss1 = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.low) ## Low smoothpar
  lines(ss1$x, ss1$y, col = 'firebrick')
  
  
  ### Plot season-specific
  par(mar = c(1.5,3,2,-.5)+1.5)
  ## Country-specific plots
  ssns = unique(splotdat$season)
  for(ss in ssns){
    valid1 = subset(splotdat, season == ss & flutype == "1")
    valid2 = subset(splotdat, season == ss & flutype == "2")
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue', percent = 70))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste(ss), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz1 = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.low)
      lines(zz1$x, zz1$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick', percent = 70))
      zz1 = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.low)
      lines(zz1$x, zz1$y, col = 'firebrick')
    }
  }
  
  
  ## Country-specific plots
  ccs = unique(cplotdat$country)
  for(cc in ccs){
    valid1 = subset(cplotdat, country == cc & flutype == "1")
    valid2 = subset(cplotdat, country == cc & flutype == "2")
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text =paste(cc), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz1 = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.low)
      lines(zz1$x, zz1$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick', percent = 70))
      zz1 = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.low)
      lines(zz1$x, zz1$y, col = 'firebrick')
    }
  }
  dev.off()
}




##################################
## Plot all countries
##################################
{
pdf(outfile4, height = 9)

par(mfrow = c(4,4))
ccs = unique(dat.002$country) # GEt all country names
## Exclude Portugal, which only reported negative cases
ccs = ccs[-7]
labs = LETTERS[1:20]; names(labs) = c(ccs)

### Plot country-specific
par(mar = c(1.5,2.5,2,0)+1.5)
## Country-specific plots
for(cc in ccs){
  valid1 = subset(dat.002, country == cc & flutype == "1")
  valid2 = subset(dat.002, country == cc & flutype == "2")
  ## Determine y limits for each plot
  mxs = c(age.tab(valid1$age)/nrow(valid1), age.tab(valid2$age)/nrow(valid2)); mxs = mxs[which(is.finite(mxs))] # Extract finite proportions and then use the maximum to set the y limits
  ymax = max(mxs)*1.05
  
  ## Plot H1N1 age dist
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, ymax), main = '', col = tns('dodgerblue'), bty = 'n')
  mtext(text = 'age', side = 1, line = 1.9, cex = .8)
  mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
  mtext(text =paste(cc), side = 3, line = .25, font = 1, cex = .9)
  mtext(labs[cc], side = 3, line = 1.2, at = 0, font = 2)
  legend(43, par("usr")[4]*.9, legend = paste(c('H1N1: n=', 'H3N2: n='), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n', xpd = NA, cex = .9)
  if(nrow(valid1) >50){# If cases of this subtype were observed in the country of interest...
  ## Add smoothing spline
  ss = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
  lines(ss$x, ss$y, col = 'dodgerblue')
  }
  if(nrow(valid2)>0){
    ## Plot H3N2 age distribution
    points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
  }
  if(nrow(valid2)>50){ # Add a smoothing spline
  ss = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
  lines(ss$x, ss$y, col = 'firebrick')
  }
}
dev.off()
}



##################################
## Plot all seasons
##################################
{
  pdf(outfile5, height = 9)
  
  par(mfrow = c(4,4))
  ssns = unique(dat.002$season) # GEt all country names
  ssns = ssns[-14] # #xclude SH.17 where no confirmed cases were observed
  ## Exclude Portugal, which only reported negative cases
  labs = LETTERS[1:20]; names(labs) = c(ssns)
  
  ### Plot country-specific
  par(mar = c(1.5,2.5,2,0)+1.5)
  ## Country-specific plots
  for(cc in ssns){
    valid1 = subset(dat.002, season == cc & flutype == "1")
    valid2 = subset(dat.002, season == cc & flutype == "2")
    ## Determine y limits for each plot
    mxs = c(age.tab(valid1$age)/nrow(valid1), age.tab(valid2$age)/nrow(valid2)); mxs = mxs[which(is.finite(mxs))] # Extract finite proportions and then use the maximum to set the y limits
    ymax = max(mxs)*1.1
    ## Plot H1N1 age dist
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, ymax), main = '', col = tns('dodgerblue'), bty = 'n')
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text =paste(cc), side = 3, line = .25, font = 1, cex = .9)
    mtext(labs[cc], side = 3, line = 1.2, at = 0, font = 2)
    legend(43, par("usr")[4]*1.05, legend = paste(c('H1N1: n=', 'H3N2: n='), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n', xpd = NA, cex = .9)
    if(nrow(valid1) >50){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      ss = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
      lines(ss$x, ss$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age distribution
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
    }
    if(nrow(valid2)>50){ # Add a smoothing spline
      ss = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
      lines(ss$x, ss$y, col = 'firebrick')
    }
  }
  dev.off()
}









###################################
## Repeat for AZ data
###################################
setwd('2017_AZ/')## Change directories to import AZ data
source('00-Inputs_multinomial.R')
setwd('../') ## Change back

## Rename seasons for plotting
raw.dat$season = gsub(pattern = '(\\d{4})(\\d{2})', replacement = '\\1-20\\2', x = raw.dat$season)

### Plot age distributions
agerange = 0:90
## Define a function to tabluate 
age.tab = function(agevec){sapply(agerange, FUN = function(aa) sum(agevec == aa))}


## Supplementary figure
##  Plot data from all seasons in which 10 or more cases of a given subtype were observed
##  Plot splines only if 50 or more cases were observed
tb = (table(raw.dat$season, raw.dat$subtype)) # Get counts for each country, H1N1 = 1, H3N2 = 2
include = rownames(tb)[which(rowSums(tb)>0)]
include = include[-which(include %in% c('2008-2009', '2009-2010'))]

## Extract data from countries to be included
azplotdat = subset(raw.dat, season %in% include & age %in% agerange)
ssns = unique(azplotdat$season)

{
  pdf(outfile6, width = 6.5, height = 10)
  ## Set smoothing paramter
  par(mfrow = c(5,3))
  par(mar = c(2,1.5,2,1)+1.5)
  labs = LETTERS[1:15]; names(labs) = ssns
  for(ss in ssns){
    valid1 = subset(azplotdat, season == ss & subtype == "H1")
    valid2 = subset(azplotdat, season == ss & subtype == "H3")
    ## Determine y limits for each plot
    if(nrow(valid1) > 10 & nrow(valid2)>10){
    mxs = c(age.tab(valid1$age)/nrow(valid1), age.tab(valid2$age)/nrow(valid2))
    }else if (nrow(valid1) > 10){
      mxs = age.tab(valid1$age)/nrow(valid1)
    }else if (nrow(valid2) > 10){
      mxs = age.tab(valid2$age)/nrow(valid2)
    }
    mxs = mxs[which(is.finite(mxs))] # Extract finite proportions and then use the maximum to set the y limits
    ymax = max(mxs)*1.1
    ## If less than 10 cases observd, plot white
    c.in = ifelse(test = nrow(valid1)>10, tns('dodgerblue'), 'white')
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(0, 90), ylim = c(0, ymax), main = '', col = c.in)
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste('AZ,', ss), side = 3, line = .25, font = 1, cex = .9)
    mtext(labs[ss], side = 3, line = 1.2, at = -10, font = 2)
    legend(20, par('usr')[4]*.9, legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >50){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
      lines(zz$x, zz$y, col = 'dodgerblue')
    }
    ## Plot H3N2 age spline
    c.in = ifelse(test = nrow(valid2)>10, tns('firebrick'), tns('white', percent = 100))
    points(agerange, age.tab(valid2$age)/nrow(valid2), col = c.in)
    if(nrow(valid2)>50){
      zz = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
      lines(zz$x, zz$y, col = 'firebrick')
    }
  }
  dev.off()
}





## Alternate version of the main text figure, except
## Plot with all ages (0-90)
## Extract data from seasons wtih more than 50 observations of each subtype
{
  pdf(outfile7, width = 6.5, height = 4.5)
  ## Set smoothing paramter
  #layout(matrix(c(1,8,2, 9,3,4, 5,6,7), byrow = T, nrow = 3))
  par(mfrow = c(2,3),mar = c(1,3,2.5,0))
  ## Overall
  valid1 = subset(azplotdat,subtype == "H1")
  valid2 = subset(azplotdat,subtype == "H3")
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(0, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'), xpd = NA)
  mtext(text = 'age', side = 1, line = 1.9, cex = .8)
  mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
  mtext(text = 'Overall', side = 3, line = .25, font = 1)
  mtext(labs[1], side = 3, line = 1.2, at = -10, font = 2)
  #text(x = 75, y = .08, paste( 'H1N1, n = ', nrow(valid1), '\nH3N2, n = ', nrow(valid2)), cex = .9)
  legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
  ## add smoothed density
  ss = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
  lines(ss$x, ss$y, col = 'dodgerblue')
  points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
  ss = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
  lines(ss$x, ss$y, col = 'firebrick')
  
  ## Extract focal seasons
  tb = (table(raw.dat$season, raw.dat$subtype)) # Get counts for each country, H1N1 = 1, H3N2 = 2
  include = rownames(tb)[which(tb[,1]>=50 & tb[,2]>=50)]
  include = include[-which(include %in% c('2008-2009', '2009-2010'))] # exclude pandemic seasons
  ## Extract data from countries to be included
  azplotdat = subset(raw.dat, season %in% include & age %in% agerange)
  ssns = unique(azplotdat$season)
  labs = LETTERS[1:10]; names(labs) = c(NA, ssns)
  
  
  par(mar = c(2,2,2,0)+1.5)
  for(ss in ssns){
    valid1 = subset(azplotdat, season == ss & subtype == "H1")
    valid2 = subset(azplotdat, season == ss & subtype == "H3")
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(0, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste('AZ,', ss), side = 3, line = .25, font = 1, cex = .9)
    mtext(labs[ss], side = 3, line = 1.2, at = -10, font = 2)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
      lines(zz$x, zz$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
      zz = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
      lines(zz$x, zz$y, col = 'firebrick')
    }
  }
  dev.off()
}






########## repeat with a higher and lower smoothing paramter (sparlow = .6, sparhigh = 1)
{
  pdf(outfile9, width = 6.5, height = 8.5)
  ## Set smoothing paramter
  par(mfrow = c(4,3))
  ## Overall
  valid1 = subset(azplotdat,subtype == "H1")
  valid2 = subset(azplotdat,subtype == "H3")
 
 ################ Plot season-specific observations with a low smoothing par 
  par(mar = c(2,1.5,2,1)+1.5)
  labs = LETTERS[1:6]; names(labs) = ssns
  for(ss in ssns){
    valid1 = subset(azplotdat, season == ss & subtype == "H1")
    valid2 = subset(azplotdat, season == ss & subtype == "H3")
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(0, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste('AZ,', ss), side = 3, line = .25, font = 1, cex = .9)
    mtext(labs[ss], side = 3, line = 1.2, at = -10, font = 2)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.low)
      lines(zz$x, zz$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
      zz = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.low)
      lines(zz$x, zz$y, col = 'firebrick')
    }
  }
  plot.new()
  ## Repeat with a high smoothpar value
  labs = LETTERS[7:12]; names(labs) = ssns
  for(ss in ssns){
    valid1 = subset(azplotdat, season == ss & subtype == "H1")
    valid2 = subset(azplotdat, season == ss & subtype == "H3")
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(0, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste('AZ,', ss), side = 3, line = .25, font = 1, cex = .9)
    mtext(labs[ss], side = 3, line = 1.2, at = -10, font = 2)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.high)
      lines(zz$x, zz$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
      zz = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.high)
      lines(zz$x, zz$y, col = 'firebrick')
    }
  }
  dev.off()
}





##### repeat wtih ages 18:90, which matches the INSIGHT data
agerange = 18:90
## Extract data from seasons wtih more than 50 observations of each subtype
tb = (table(raw.dat$season, raw.dat$subtype)) # Get counts for each country, H1N1 = 1, H3N2 = 2
include = rownames(tb)[which(tb[,1]>=50 & tb[,2]>=50)]
include = include[-2] # Exclude pandemic seasons
## Extract data from countries to be included
azplotdat = subset(raw.dat, season %in% include & age %in% agerange)
ssns = unique(azplotdat$season)

{
  pdf(outfile8, width = 6.5, height = 6)
  ## Set smoothing paramter
  #layout(matrix(c(1,8,2, 9,3,4, 5,6,7), byrow = T, nrow = 3))
  par(mfrow = c(2,3), mar = c(1.5,3,2.5,0))
  ## Overall
  valid1 = subset(azplotdat,subtype == "H1")
  valid2 = subset(azplotdat,subtype == "H3")
  plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'), xpd = NA)
  mtext(text = 'age', side = 1, line = 1.9, cex = .8)
  mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
  mtext(text = 'Overall', side = 3, line = .25, font = 2)
  #text(x = 75, y = .08, paste( 'H1N1, n = ', nrow(valid1), '\nH3N2, n = ', nrow(valid2)), cex = .9)
  legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
  ## add smoothed density
  ss = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
  lines(ss$x, ss$y, col = 'dodgerblue')
  points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
  ss = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
  lines(ss$x, ss$y, col = 'firebrick')
  
  par(mar = c(2,2,2,0)+1.5)
  for(ss in ssns){
    valid1 = subset(azplotdat, season == ss & subtype == "H1")
    valid2 = subset(azplotdat, season == ss & subtype == "H3")
    plot(agerange, age.tab(valid1$age)/nrow(valid1), xlab = '', ylab = '', xlim = c(18, 90), ylim = c(0, .085), main = '', col = tns('dodgerblue'))
    mtext(text = 'age', side = 1, line = 1.9, cex = .8)
    mtext(text = 'fraction', side = 2, line = 1.9, cex = .8)
    mtext(text = paste(ss), side = 3, line = .25, font = 2, cex = .9)
    legend('topright', legend = paste(c('H1N1: n = ', 'H3N2: n = '), c(nrow(valid1), nrow(valid2)), sep = ''), col = c('dodgerblue', 'firebrick'), pch = 15, bty = 'n')
    if(nrow(valid1) >0){# If cases of this subtype were observed in the country of interest...
      ## Add smoothing spline
      zz = smooth.spline(agerange, age.tab(valid1$age)/nrow(valid1), spar = spar.in)
      lines(zz$x, zz$y, col = 'dodgerblue')
    }
    if(nrow(valid2)>0){
      ## Plot H3N2 age spline
      points(agerange, age.tab(valid2$age)/nrow(valid2), col = tns('firebrick'))
      zz = smooth.spline(agerange, age.tab(valid2$age)/nrow(valid2), spar = spar.in)
      lines(zz$x, zz$y, col = 'firebrick')
    }
  }
  dev.off()
}



