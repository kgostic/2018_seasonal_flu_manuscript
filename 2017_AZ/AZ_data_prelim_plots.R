#### AZ Dept. of public health H1, H3 surveillance ####


# Set directory
setwd('~/Dropbox/R/2017_seasonal_flu/')
# Clear memory
rm(list = ls())
library(ggplot2)
library(reshape)
source('multiplot.R')
# Load data
filenames = paste('AZ_flu_surveillance_', c('92', '93', '94', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14'), '-',
                  c('93', '94', '95', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'), sep = '')
#dat.92.93 = read.csv('AZ_flu_surveillance_92-93.csv', stringsAsFactors = FALSE) # Import as a separate variable, since only contains one data point
all.data = read.csv('AZ_flu_surveillance_93-94.csv', stringsAsFactors = FALSE)
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_94-95.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_02-03.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_03-04.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_04-05.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_05-06.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_06-07.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_07-08.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_08-09.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_09-10.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_10-11.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_11-12.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_12-13.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_13-14.csv', stringsAsFactors = FALSE))
all.data = rbind(all.data, read.csv('AZ_flu_surveillance_14-15.csv', stringsAsFactors = FALSE))

# drop 2046 birth year typo
all.data = all.data[-which(all.data$BIRTHYEAR == 2046), ]


# Add age to all.data
current.year = as.numeric(gsub("(\\d\\d\\d\\d)(\\d\\d)", "\\1", as.character(all.data$SEASON)))
all.data$age = current.year - all.data$BIRTHYEAR +1
all.data$year = current.year

## Break into H1 and H3 subsets
H3 = subset(all.data, FLU_SUBTYPE == 'H3')
H1 = subset(all.data, FLU_SUBTYPE == 'H1')


## Save data
write.csv(H1, file = 'AZ_H1.csv')
write.csv(H3, file = 'AZ_H3.csv')



#Sort by season and birth year
H1.table = table(H1[,c(1,2)])
H3.table = table(H3[,c(1,2)])



#Plot raw case counts
plot1  = ggplot(H1, aes(x = BIRTHYEAR)) + geom_histogram(aes(fill = as.factor(SEASON)), binwidth = 1) + labs(title = 'H1N1', x = 'Birth year', y = 'Cases') + scale_x_continuous(breaks=seq(1918, 2018, 10)) + theme(axis.text=element_text(size=12))

plot2  = ggplot(H3, aes(x = BIRTHYEAR)) + geom_histogram(aes(fill = as.factor(SEASON)), binwidth = 1) + labs(title = 'H3N2', x = 'Birth year', y = 'Cases') + scale_x_continuous(breaks=seq(1918, 2018, 10)) + theme(axis.text=element_text(size=12))

pdf('RawCaseCounts.pdf')
multiplot(plot1, plot2)
dev.off()





## Base Graphics Barplots
pdf('H1raw.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
cols = cols.H1=  (c('deeppink', 'palevioletred1', 'green3', 'palegreen1', 'seagreen', 'mediumblue', 'dodgerblue', 'lightskyblue', 'navy', 'mediumpurple4', 'slateblue', 'darkmagenta'))
xx = barplot(H1.table, names.arg = colnames(H1.table), col = cols, ylab = 'Cases', main = 'H1N1', border = cols, xlab = 'Birth year')
#Add demographic curve
source('Get_US_census.R')
# Normalize by observation year
dem.curve.raw = (US.demography/rowSums(US.demography))[c(2, 3, 6:15), as.character(1918:2015)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H1.table))
use.H1.table = which(colnames(H1.table) %in% colnames(dem.curve.raw))
dem.curve.raw = dem.curve.raw[, use.dem.curve]
legend.labs = gsub("(\\d\\d\\d\\d)(\\d\\d)", "\\1-\\2", rownames(H1.table))
legend('topleft', legend = legend.labs, fill = cols, ncol = 2, bty = 'n')
H1.table = H1.table[, use.H1.table]
dem.curve = colSums(dem.curve.raw*rowSums(H1.table))
lines(xx[use.H1.table], (dem.curve), lwd = 3)
dev.off()


## Base Graphics Barplots
pdf('H3raw.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
cols = cols.H3 = (c('red', 'deeppink', 'palevioletred1', 'orange', 'gold', 'green3', 'palegreen1', 'seagreen', 'mediumblue', 'dodgerblue', 'lightskyblue', 'navy', 'mediumpurple4', 'slateblue', 'darkmagenta'))
xx = barplot(H3.table, names.arg = colnames(H3.table), col = cols, ylab = 'Cases', main = 'H3N2', border = cols, xlab = 'Brith year')
leg.labs = gsub("(\\d\\d\\d\\d)(\\d\\d)", "\\1-\\2", rownames(H3.table))
legend('topleft', legend = leg.labs, fill = cols, ncol = 2, bty = 'n')
#Add demographic curve
source('Get_US_census.R')
# Normalize by observation year
dem.curve.raw = (US.demography/rowSums(US.demography))[, as.character(1918:2015)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H3.table))
use.H3.table = which(colnames(H3.table) %in% colnames(dem.curve.raw))
dem.curve.raw = dem.curve.raw[, use.dem.curve]
H3.table = H3.table[, use.H3.table]
dem.curve = colSums(dem.curve.raw*rowSums(H3.table))
lines(xx[use.H3.table], (dem.curve), lwd = 3)
dev.off()



## Get imprinting  data
source('Get_US_vax_weights.R')
H1.mismatched.weights = weights.master.3 + weights.master.naiive
H1.mismatched.weights = H1.mismatched.weights[which(rownames(H1.mismatched.weights) %in% c('1995USA', '2003USA', '2006USA', '2007USA', '2008USA', '2009USA', '2010USA', '2011USA', '2012USA', '2013USA', '2014USA', '2015USA')), as.character(1918:2014)]
H3.mismatched.weights = weights.master.1+weights.master.2+weights.master.naiive
H3.mismatched.weights = H3.mismatched.weights[which(rownames(H3.mismatched.weights) %in% c('1994USA', '1995USA', '2003USA', '2004USA', '2005USA', '2006USA', '2007USA', '2008USA', '2009USA', '2010USA', '2011USA', '2012USA', '2013USA', '2014USA', '2015USA')), as.character(1918:2015)]




## Normalized
pdf('H1norm.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
source('Get_US_census.R')
dem.curve.raw = (US.demography/rowSums(US.demography))[c(2, 3, 6:15), as.character(1918:2014)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H1.table))
use.H1.table = which(colnames(H1.table) %in% colnames(dem.curve.raw))
exp.H1 = (dem.curve.raw[,use.dem.curve]*rowSums(H1.table)) #This was dem.curve last time
obs.H1 = H1.table[,use.H1.table]
xx = barplot((obs.H1-exp.H1), names.arg = colnames(obs.H1), col = cols.H1, ylab = 'Excess Cases', main = 'H1N1', border = cols.H1, xlab = 'Birth year', ylim = c(-120, 300))
legend('topleft', legend = c(legend.labs), fill = c(cols.H1), ncol = 2, bty = 'n')
imprinting.H1 = H1.mismatched.weights*dem.curve.raw/rowSums(H1.mismatched.weights*dem.curve.raw)
par(mar = c(2, 7, 3, 3), mfrow =c (3,1))
plot(1918:2014, H1.mismatched.weights[12,], lwd = 2, type = 'l', ylab = '', xlab = 'Birth year')
plot(1918:2015, H3.mismatched.weights[15, ], lwd = 2, type = 'l', ylab = '', xlab = 'Birth year')
dev.off()


pdf('H1norm_withimp.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
source('Get_US_census.R')
dem.curve.raw = (US.demography/rowSums(US.demography))[c(2, 3, 6:15), as.character(1918:2014)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H1.table))
use.H1.table = which(colnames(H1.table) %in% colnames(dem.curve.raw))
exp.H1 = (dem.curve.raw[,use.dem.curve]*rowSums(H1.table)) #This was dem.curve last time
obs.H1 = H1.table[,use.H1.table]
xx = barplot((obs.H1-exp.H1), names.arg = colnames(obs.H1), col = cols.H1, ylab = 'Excess Cases', main = 'H1N1', border = cols.H1, xlab = 'Birth year', ylim = c(-200, 450))
imprinting.H1 = H1.mismatched.weights*dem.curve.raw/rowSums(H1.mismatched.weights*dem.curve.raw)*rowSums(H1.table)
lines(xx, colSums(imprinting.H1-exp.H1), lwd = 2)
legend('topleft', legend = c(legend.labs), fill = c(cols.H1), ncol = 3, bty = 'n')
dev.off()

## Base Graphics Barplots
pdf('H1raw_withimp.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
cols = cols.H1=  (c('deeppink', 'palevioletred1', 'green3', 'palegreen1', 'seagreen', 'mediumblue', 'dodgerblue', 'lightskyblue', 'navy', 'mediumpurple4', 'slateblue', 'darkmagenta'))
xx = barplot(H1.table, names.arg = colnames(H1.table), col = cols, ylab = 'Cases', main = 'H1N1', border = cols, xlab = 'Birth year')
#Add demographic curve
source('Get_US_census.R')
# Normalize by observation year
dem.curve.raw = (US.demography/rowSums(US.demography))[c(2, 3, 6:15), as.character(1918:2015)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H1.table))
use.H1.table = which(colnames(H1.table) %in% colnames(dem.curve.raw))
dem.curve.raw = dem.curve.raw[, use.dem.curve]
legend.labs = gsub("(\\d\\d\\d\\d)(\\d\\d)", "\\1-\\2", rownames(H1.table))
legend('topleft', legend = legend.labs, fill = cols, ncol = 1, bty = 'n')
H1.table = H1.table[, use.H1.table]
dem.curve = colSums(dem.curve.raw*rowSums(H1.table))
lines(xx[use.H1.table], colSums(imprinting.H1), lwd = 3)
dev.off()


## Base Graphics Barplots
pdf('H3norm.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
source('Get_US_census.R')
dem.curve.raw = (US.demography/rowSums(US.demography))[, as.character(1918:2015)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H3.table))
use.H3.table = which(colnames(H3.table) %in% colnames(dem.curve.raw))
exp.H3 = (dem.curve.raw[, use.dem.curve]*rowSums(H3.table))
obs.H3 = H3.table[,use.H3.table]
xx = barplot((obs.H3-exp.H3), names.arg = colnames(obs.H3), col = cols.H3, ylab = 'Excess Cases', main = 'H3N2', border = cols, xlab = 'Birth year', ylim = c(-50, 200))
legend('topleft', legend = leg.labs, fill = cols, ncol = 2, bty = 'n')
dev.off()



## Base Graphics Barplots
pdf('H3norm_imp.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
source('Get_US_census.R')
dem.curve.raw = (US.demography/rowSums(US.demography))[, as.character(1918:2015)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H3.table))
use.H3.table = which(colnames(H3.table) %in% colnames(dem.curve.raw))
exp.H3 = (dem.curve.raw[, use.dem.curve]*rowSums(H3.table))
obs.H3 = H3.table[,use.H3.table]
xx = barplot((obs.H3-exp.H3), names.arg = colnames(obs.H3), col = cols.H3, ylab = 'Excess Cases', main = 'H3N2', border = cols, xlab = 'Birth year', ylim = c(-100, 200))
imprinting.H3 = H3.mismatched.weights*dem.curve.raw/rowSums(H3.mismatched.weights*dem.curve.raw)*rowSums(H3.table)
lines(xx, colSums(imprinting.H3-exp.H3), lwd = 2)
legend('topleft', legend = c(leg.labs), fill = c(cols.H3), ncol = 3, bty = 'n')
dev.off()




pdf('H1raw_HAimp.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
cols = rev(c('red', 'deeppink', 'palevioletred1', 'orange', 'gold', 'green3', 'mediumblue', 'dodgerblue', 'lightskyblue', 'navy', 'mediumpurple4'))
xx = barplot(H1.table, names.arg = colnames(H1.table), col = cols, ylab = 'Cases', main = 'H1N1', border = cols, xlab = 'Birth year', ylim = c(0, 600))
legend('topleft', legend = rownames(H1.table), fill = cols, ncol = 2, cex  = .9)
#Add demographic curve
source('../ImmuneAgeStructure_ch2/Get_US_census.R')
# Normalize by observation year
H1.wts = colMeans(H1.mismatched.weights)*300
use.H1.wts = which(names(H1.wts) %in% colnames(H1.table))
use.H1.table = which(colnames(H1.table) %in% names(H1.wts))
lines(xx[use.H1.table], rev(H1.wts[use.H1.wts]), lwd = 3)
lines(xx[use.H1.table], rev(H1.wts[use.H1.wts])-3, lwd = 2, col = 'white')
dev.off()

## Base Graphics Barplots
pdf('H3raw_HAimp.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
cols = (c('red', 'deeppink', 'palevioletred1', 'orange', 'gold', 'green3', 'palegreen1', 'seagreen', 'mediumblue', 'dodgerblue', 'lightskyblue', 'navy', 'mediumpurple4', 'slateblue', 'darkmagenta'))
xx = barplot(H3.table, names.arg = colnames(H3.table), col = cols, ylab = 'Cases', main = 'H3N2', border = cols, xlab = 'Brith year', ylim = c(0, 350))
legend('topleft', legend = rownames(H3.table), fill = cols, ncol = 2, cex = .8)
#Add demographic curve
source('~/Dropbox/R/2017_seasonal_flu/Get_US_census.R')
H3.wts = colMeans(H3.mismatched.weights)*150
use.H3.wts = which(names(H3.wts) %in% colnames(H3.table))
use.H3.table = which(colnames(H3.table) %in% names(H3.wts))
lines(xx[use.H3.table], rev(H3.wts[use.H3.wts]), lwd = 3)
lines(xx[use.H3.table], rev(H3.wts[use.H3.wts])-3, lwd = 2, col = 'white')
# Normalize by observation year
dev.off()


## Base Graphics Barplots
pdf('H3raw_withimp.pdf')
par(cex = 1.4, cex.axis = 1.4, cex.lab = 1.4)
cols = cols.H3 = (c('red', 'deeppink', 'palevioletred1', 'orange', 'gold', 'green3', 'palegreen1', 'seagreen', 'mediumblue', 'dodgerblue', 'lightskyblue', 'navy', 'mediumpurple4', 'slateblue', 'darkmagenta'))
xx = barplot(H3.table, names.arg = colnames(H3.table), col = cols, ylab = 'Cases', main = 'H3N2', border = cols, xlab = 'Brith year')
leg.labs = gsub("(\\d\\d\\d\\d)(\\d\\d)", "\\1-\\2", rownames(H3.table))
legend('topleft', legend = leg.labs, fill = cols, ncol = 2, bty = 'n')
#Add demographic curve
source('Get_US_census.R')
# Normalize by observation year
dem.curve.raw = (US.demography/rowSums(US.demography))[, as.character(1918:2015)]
use.dem.curve = which(colnames(dem.curve.raw) %in% colnames(H3.table))
use.H3.table = which(colnames(H3.table) %in% colnames(dem.curve.raw))
dem.curve.raw = dem.curve.raw[, use.dem.curve]
H3.table = H3.table[, use.H3.table]
dem.curve = colSums(dem.curve.raw*rowSums(H3.table))
lines(xx[use.H3.table], colSums(imprinting.H3), lwd = 3)
dev.off()





#List seasons
seasons = c(199394, 199495, 200203, 200304, 200405, 200506, 200607, 200708, 200809, 200910, 201011, 201112, 201213, 201314, 201415)
#Initialize data frames to store the fraction of cases in each birth year by season
norm.H1 = norm.H3 = data.frame(BIRTHYEAR = rep(1903:2015, length(seasons)), case.fraction = NA, season = rep(seasons, each = length(1903:2015)))

# Fill in by season
for(ii in 1:length(seasons)){
  if(seasons[ii] %in% H1$SEASON){ # Skip if no H1 cases observed in season ii
  season.dat.H1 = subset(H1, SEASON == seasons[ii]) # If cases observed, extract relevant data
  if(nrow(season.dat.H1) > 30){ # Fill in normalized data only if >30 observations, else leave NA
    # Count up the number of entries in each birth year
    norm.H1[norm.H1$season == seasons[ii], 'case.fraction'] = sapply(1903:2015, function(x) sum(season.dat.H1$BIRTHYEAR == x))/length(season.dat.H1$BIRTHYEAR)
  }
  }
  
  #repeat for H3
  if(seasons[ii] %in% H3$SEASON){ # Skip if no H3 cases observed in season ii
    season.dat.H3 = subset(H3, SEASON == seasons[ii]) # If cases observed, extract relevant data
    if(nrow(season.dat.H3) > 30){ # Fill in normalized data only if >30 observations, else leave NA
      # Count up the number of entries in each birth year
      norm.H3[norm.H3$season == seasons[ii], 'case.fraction'] = sapply(1903:2015, function(x) sum(season.dat.H3$BIRTHYEAR == x))/length(season.dat.H3$BIRTHYEAR)
    }
  }
  
}


    
  
  
#Plot normalized case counts
plot1.norm  = ggplot(norm.H1, aes(x = (BIRTHYEAR), y = case.fraction)) + geom_histogram(stat = "identity", aes(fill = as.factor(season)), binwidth = 1) + labs(title = 'H1N1 normalized', x = 'Birth year', y = 'Cases') + scale_x_continuous(breaks=seq(1918, 2018, 10)) + theme(axis.text=element_text(size=12))


plot2.norm = ggplot(norm.H3, aes(x = (BIRTHYEAR), y = case.fraction)) + geom_histogram(stat = "identity", aes(fill = as.factor(season)), binwidth = 1) + labs(title = 'H3N2 normalized', x = 'Birth year', y = 'Cases') +  scale_x_continuous(breaks=seq(1918, 2018, 10)) + theme(axis.text=element_text(size=12))

pdf('Normalized_Case_Counts.pdf')
multiplot(plot1.norm, plot2.norm)
dev.off()

pdf('Raw_and_Norm_CaseCounts.pdf')
multiplot(plot1, plot2, plot1.norm, plot2.norm, cols = 2)
dev.off()


# Plot case observations by season
par(mfrow = c(3,3), mar = c(2, 2, 1, 2))
for(ii in 1:length(unique(norm.H1$season))){
  if( sum(is.na(subset(norm.H1, season == unique(season)[ii], select = case.fraction))) < 10){
  barplot(height =  subset(norm.H1, season == unique(season)[ii], select = case.fraction)$case.fraction, names.arg = 1903:2015, main = sprintf('H1N1 %s', seasons[ii]))
  }
}

# Plot case observations by season
par(mfrow = c(4,3), mar = c(2, 2, 1, 2))
for(ii in 1:length(unique(norm.H3$season))){
  if( sum(is.na(subset(norm.H3, season == unique(season)[ii], select = case.fraction))) < 10){
    barplot(height =  subset(norm.H3, season == unique(season)[ii], select = case.fraction)$case.fraction, names.arg = 1903:2015, main = sprintf('H3N2, %s', seasons[ii]))
  }
}



## First remove seasons with < 30 samples
norm.H1.clean = norm.H1[-which(is.na(norm.H1$case.fraction)), ]
norm.H3.clean = norm.H3[-which(is.na(norm.H3$case.fraction)), ]

## plot normalized fraction of cases for H1
pp1 = ggplot(norm.H1.clean, aes(x = BIRTHYEAR, y = case.fraction))
pp1 + geom_point(aes(color = as.factor(season))) + ggtitle('H1N1')

## Get imprinting  data
source('Get_US_vax_weights.R')
H1.mismatched.weights = weights.master.3 + weights.master.naiive
H1.mismatched.weights = H1.mismatched.weights[-which(rownames(H1.mismatched.weights) %in% c('1995USA', '1996USA', '1997USA', '1998USA', '1999USA', '2000USA', '2001USA','2015USA')), ]
H3.mismatched.weights = weights.master.1+weights.master.2+weights.master.naiive
H3.mismatched.weights = H3.mismatched.weights[-which(rownames(H3.mismatched.weights) %in% c('1995USA', '1996USA', '1997USA', '1998USA', '1999USA', '2000USA', '2001USA', '2015USA')), ]

# Replace rownames to match seasons
rownames(H1.mismatched.weights) = rownames(H3.mismatched.weights) = seasons

# Melt the data frame for ggplot
H1.mismatched.weights = melt(H1.mismatched.weights); colnames(H1.mismatched.weights) = c('season', 'birthyear', 'mismatch.frac')
H3.mismatched.weights = melt(H3.mismatched.weights); colnames(H3.mismatched.weights) = c('season', 'birthyear', 'mismatch.frac')

# Reorder by birth year
H1.mismatched.weights = H1.mismatched.weights[order(H1.mismatched.weights$season, H1.mismatched.weights$birthyear), ]
H3.mismatched.weights = H3.mismatched.weights[order(H3.mismatched.weights$season, H3.mismatched.weights$birthyear), ]

# Get rid of pre-1918 birth years
norm.H1.1918on = norm.H1[-which(norm.H1$BIRTHYEAR < 1918), ]
norm.H3.1918on = norm.H3[-which(norm.H3$BIRTHYEAR < 1918), ]

# Combine data frames
H1.imprinting.correlation = data.frame(season = H1.mismatched.weights$season, birthyear = H1.mismatched.weights$birthyear, observed.case.frac = norm.H1.1918on$case.fraction, mismatch.frac = H1.mismatched.weights$mismatch.frac)

H3.imprinting.correlation = data.frame(season = H3.mismatched.weights$season, birthyear = H3.mismatched.weights$birthyear, observed.case.frac = norm.H3.1918on$case.fraction, mismatch.frac = H3.mismatched.weights$mismatch.frac)


#Plot the fraction mismatched imprinting vs. observed case fraction
pp2 = ggplot(H1.imprinting.correlation, aes(x = mismatch.frac, y = observed.case.frac))
pdf('Imprinting_correlation.pdf')
pp2 + geom_jitter(aes(color = as.factor(season)), width = .5, height = .5) + geom_smooth(method = lm)

pp3 = ggplot(H3.imprinting.correlation, aes(x = mismatch.frac, y = observed.case.frac))
pp3 + geom_jitter(aes(color = as.factor(season))) + geom_smooth(method = lm)
dev.off()




cor(H1.imprinting.correlation$observed.case.frac, H1.imprinting.correlation$mismatch.frac, use = "complete.obs")

# H1 correlation 
cor.H1 = na.omit(H1.imprinting.correlation)
cor.test(cor.H1$observed.case.frac, cor.H1$mismatch.frac, method = "spearman")

# H3 correlation
cor.H3 = na.omit(H3.imprinting.correlation)
cor.test(cor.H3$observed.case.frac, cor.H3$mismatch.frac, method = "spearman")


ggplot(data = cor.H1[with(cor.H1, order(mismatch.frac)), ], aes(y = observed.case.frac, x = mismatch.frac)) + geom_jitter(aes(color = as.factor(birthyear))) + geom_smooth(method = lm) + ggtitle('H1N1 Rank Correlation')

ggplot(data = cor.H3[with(cor.H3, order(mismatch.frac)), ], aes(y = observed.case.frac, x = mismatch.frac)) + geom_jitter(aes(color = as.factor(birthyear))) + geom_smooth(method = lm) + ggtitle('H3N2 Rank Correlation')




# Obtain demographic excess cases
setwd('~/Dropbox/R/ImmuneAgeStructure_ch2/')
source('Get_US_census.R')
# Normalize by observation year
US.demography = US.demography/rowSums(US.demography)

# Replace rownames to match seasons
rownames(US.demography) = seasons

# Melt the data frame for ggplot
US.demography = melt(US.demography); colnames(US.demography) = c('season', 'birthyear', 'pop.frac')
# Reorder by birth year
US.demography = US.demography[order(US.demography$season, US.demography$birthyear), ]


# Add to the H3.imprinting.correlation and H1.imprintnig.correlation data frame
case.counts = c(0, rowSums(H1.table)[1:2], 0, 0, rowSums(H1.table)[3:10], 0, 0)
n.H1.cases = rep(case.counts, each = 98)

H1.imprinting.correlation = data.frame(H1.imprinting.correlation, pop.frac = US.demography$pop.frac, excess.cases = (H1.imprinting.correlation$observed.case.frac - US.demography$pop.frac))


#How many cases occurred per year?
n.H3.cases = rep(rowSums(H3.table), each = 98)

H3.imprinting.correlation = data.frame(H3.imprinting.correlation, pop.frac = US.demography$pop.frac, excess.cases = (H3.imprinting.correlation$observed.case.frac - US.demography$pop.frac)*n.H3.cases)


# H1 correlation with demographic normalization
case.counts = rep(rowSums(H1.table)[which(rowSums(H1.table) > 30)], each = 98)
cor.H1 = cbind(na.omit(H1.imprinting.correlation), case.counts, case.counts*na.omit(H1.imprinting.correlation)[,'excess.cases'])
colnames(cor.H1)[8] = 'excess.case.counts'

case.counts = rep(rowSums(H3.table)[which(rowSums(H3.table) > 30)], each = 98)
cor.H3 = cbind(na.omit(H3.imprinting.correlation), case.counts, case.counts*na.omit(H3.imprinting.correlation)[,'excess.cases'])
colnames(cor.H3)[8] = 'excess.case.counts'

# Birth year vs. excess cases
ggplot(cor.H1, aes(x = birthyear, y = excess.cases, color = as.factor(season))) + geom_point()












pdf('Age_excess plots.pdf')
# H1N1 with imprinting expectation
ggplot(cor.H1, aes(x = birthyear, y = excess.cases, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H1N1') + xlab('Birth Year') + ylab('Observed - Expected') + geom_line(data = subset(cor.H1, season == 201314), aes(x = birthyear, y = (mismatch.frac*0.5)-0.25))

# H3N2 with imprinting expectation
ggplot(cor.H3, aes(x = birthyear, y = excess.cases, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H3N2') + xlab('Birth Year') + ylab('Observed - Expected') + 
  geom_line(data = subset(cor.H3, season == 201415), aes(x = birthyear, y = (mismatch.frac*0.25)-0.125))

# H1N1 plain
# #fractions
# ggplot(cor.H1, aes(x = birthyear, y = excess.cases, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H1N1') + xlab('Birth Year') + ylab('Observed - Expected') 

#counts
ggplot(cor.H1, aes(x = birthyear, y = excess.case.counts, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H1N1') + xlab('Birth Year') + ylab('Observed - Expected') 
# H3N2 plain
# #fractions
# ggplot(cor.H3, aes(x = birthyear, y = excess.cases, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H3N2') + xlab('Birth Year') + ylab('Observed - Expected')

#counts
ggplot(cor.H3, aes(x = birthyear, y = excess.case.counts, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H3N2') + xlab('Birth Year') + ylab('Observed - Expected')
dev.off()



pdf('Age_excess plots2.pdf', height = 3.5)

#counts
ggplot(cor.H1, aes(x = birthyear, y = excess.case.counts, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H1N1') + xlab('Birth Year') + ylab('Excess cases') 
# H3N2 plain
# #fractions
# ggplot(cor.H3, aes(x = birthyear, y = excess.cases, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H3N2') + xlab('Birth Year') + ylab('Observed - Expected')

#counts
ggplot(cor.H3, aes(x = birthyear, y = excess.case.counts, fill = as.factor(season))) + geom_bar(stat = 'identity', bindidth = 1) + ggtitle('H3N2') + xlab('Birth Year') + ylab('Excess cases')
dev.off()


# This doesn't really meet the heteroskedacitiy assumption
ggplot(cor.H3, aes(x = birthyear, y = excess.cases, color = as.factor(season))) + geom_point()
# This does though

# H1 rank correlation of excess cases with mismatch fraction
cor.test(cor.H1$excess.cases, cor.H1$mismatch.frac, method = "spearman")

# H3 rank correlation of excess cases with mismatch fraction
cor.H3 = na.omit(H3.imprinting.correlation)
cor.test(cor.H3$excess.cases, cor.H3$mismatch.frac, method = "spearman")



ggplot(data = cor.H1[with(cor.H1, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = as.factor(birthyear))) + ggtitle('H1N1') + geom_smooth(method = lm)


ggplot(data = cor.H3[with(cor.H3, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = as.factor(birthyear))) + ggtitle('H3N2') + geom_smooth(method = lm)


# Try 5 year aga bins
birth.year.bins = matrix(c(NA, NA, 1918:2015), ncol = 5, byrow = T)


ss = subset(cor.H1, season == unique(cor.H1$season)[1])
apply(birth.year.bins, MARGIN = 1, function(x) colSums(ss[ss$birthyear %in% x, 3:6]))

#Set up 5-yr-bin matrix
cor.H1.5yrbins = as.data.frame(matrix(NA, nrow = 20*length(unique(cor.H1$season)), ncol = 6, dimnames = list(NULL, c('season', 'birthyear', 'observed.case.frac', 'mismatch.frac', 'pop.frac', 'excess.cases'))))
cor.H1.5yrbins$season = rep(unique(cor.H1$season), each = 20)
cor.H1.5yrbins$birthyear = rep( paste(c(1918, birth.year.bins[-1,1]), birth.year.bins[,5], sep = '-'), times = length(unique(cor.H1$season)))
for(ii in 1:length(unique(cor.H1$season))){
  ss = subset(cor.H1, season == unique(cor.H1$season)[ii])
  output = apply(birth.year.bins, MARGIN = 1, function(x) colSums(ss[ss$birthyear %in% x, 3:6]))
  cor.H1.5yrbins[which(cor.H1.5yrbins$season == unique(cor.H1$season)[ii]), 3:6] = t(output)
}

cor.H3.5yrbins = as.data.frame(matrix(NA, nrow = 20*length(unique(cor.H3$season)), ncol = 6, dimnames = list(NULL, c('season', 'birthyear', 'observed.case.frac', 'mismatch.frac', 'pop.frac', 'excess.cases'))))
cor.H3.5yrbins$season = rep(unique(cor.H3$season), each = 20)
cor.H3.5yrbins$birthyear = rep( paste(c(1918, birth.year.bins[-1,1]), birth.year.bins[,5], sep = '-'), times = length(unique(cor.H3$season)) )

for(ii in 1:length(unique(cor.H3$season))){
  ss = subset(cor.H3, season == unique(cor.H3$season)[ii])
  output = apply(birth.year.bins, MARGIN = 1, function(x) colSums(ss[ss$birthyear %in% x, 3:6]))
  cor.H3.5yrbins[which(cor.H3.5yrbins$season == unique(cor.H3$season)[ii]), 3:6] = t(output)
}



# H1 correlation with demographic normalization
cor.H1 = na.omit(cor.H1.5yrbins)
cor.test(cor.H1.5yrbins$excess.cases, cor.H1.5yrbins$mismatch.frac, method = "spearman")
ggplot(data = cor.H1.5yrbins[with(cor.H1.5yrbins, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = birthyear)) + ggtitle('H1N1')



H1.mismatch = colSums(weights.master.3+weights.master.naiive)/dim(weights.master.3)[1]

ggplot(cor.H1.5yrbins[-which(cor.H1.5yrbins$excess.cases == 0), ], aes(x = birthyear, y = excess.cases)) + geom_bar(stat = 'identity', aes(fill = as.factor(season))) + ggtitle('H1N1') 

yy = ggplot(data = data.frame(years = 2015:1918, weights = H1.mismatch), aes(x = years, y = weights)) + geom_line()




# H3 correlation with demographic normalization
cor.H3 = na.omit(cor.H3.5yrbins)
cor.test(cor.H3.5yrbins$excess.cases, cor.H3.5yrbins$mismatch.frac, method = "spearman")
ggplot(data = cor.H3.5yrbins[with(cor.H3.5yrbins, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = birthyear)) + ggtitle('H3N2')


ggplot(cor.H1.5yrbins, aes(x = birthyear, y = excess.cases, color = as.factor(season))) + geom_point() + geom_smooth(aes(group = as.factor(season)), se = FALSE) + ggtitle('H1N1 5 year age bins')
# This doesn't really meet the heteroskedacitiy assumption

ggplot(cor.H3, aes(x = birthyear, y = excess.cases, color = as.factor(season))) + geom_point()
# This does though


# Pull out cluster jump and non-cluster jump years
H1.jump = subset(cor.H1.5yrbins, season %in% c('200708', '200809', '201011'))
H1.no.jump = subset(cor.H1.5yrbins, !season %in% c('200708', '200809', '201011'))

H3.jump = subset(cor.H3.5yrbins, !season %in% c('199394', '200304', '200708', '200910', '201112', '201213', '201415'))
H3.no.jump = subset(cor.H3.5yrbins, season %in% c('199394','200304', '200708', '200910', '201112', '201213', '201415'))

# H3 jump
cor.test(H3.jump$excess.cases, H3.jump$mismatch.frac, method = "spearman")
cor.test(H3.no.jump$excess.cases, H3.no.jump$mismatch.frac, method = "spearman")
# Plot jump years
ggplot(data = H3.jump[with(H3.jump, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = birthyear)) + ggtitle('H3N2 cluster jumps')
# Plot no jump years
ggplot(data = H3.no.jump[with(H3.no.jump, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = birthyear)) + ggtitle('H3N2 cluster jumps')

# H1 jump
cor.test(H1.jump$excess.cases, H1.jump$mismatch.frac, method = "spearman")
cor.test(H1.no.jump$excess.cases, H1.no.jump$mismatch.frac, method = "spearman")
# Plot jump years
ggplot(data = H1.jump[with(H1.jump, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = birthyear)) + ggtitle('H1N1 cluster jumps')
# Plot no jump years
ggplot(data = H1.no.jump[with(H1.no.jump, order(mismatch.frac)), ], aes(y = excess.cases, x = mismatch.frac)) + geom_point(aes(color = birthyear)) + ggtitle('H1N1 cluster jumps')



# Find correlation for every year of circulation after 2010
spearmans_rho = numeric(length(2010:2014))
pandemic.season = c('200910', '201011', '201112', '201213', '201314')
plot(pandemic.season, sapply(pandemic.season, function(x) cor.test(~ excess.cases + mismatch.frac, data = cor.H1.5yrbins, subset = season == x, method = 'spearman')$estimate))
sapply(pandemic.season, function(x) cor.test(~ excess.cases + mismatch.frac, data = cor.H1.5yrbins, subset = season == x, method = 'spearman')$p.value)




# OTHER FACTORS -> AGE, VACCINATION, pre-existing immunity!




##### ________________________________________
####     REGRESS AGE VS. SEASON
####      BIRTH YEAR VS. SEASON
##### ________________________________________




# Plot age vs. season
par(mfrow = c(2, 1))
H1 = subset(all.data, FLU_SUBTYPE == 'H1')
plot((H1$SEASON), H1$age, ylab = 'age', xlab = 'Season', type = 'p', main = 'H1N1')
age.lm = summary(lm(age ~ SEASON, data = H1)); age.lm
lines(H1$SEASON, age.lm$coefficients[1]+age.lm$coefficients[2]*H1$SEASON)
cor.test(H1$SEASON, H1$age, method = 'spearman')

plot((H1$SEASON), H1$BIRTHYEAR, col = 'blue', ylab = 'Birth year', xlab = 'Season')
BY.lm = summary(lm(BIRTHYEAR ~ SEASON, data = H1)); BY.lm
lines(H1$SEASON, BY.lm$coefficients[1]+BY.lm$coefficients[2]*H1$SEASON)
cor.test(H1$SEASON, H1$BIRTHYEAR, method = 'spearman')


par(mfrow = c(2, 1))
H3 = subset(all.data, FLU_SUBTYPE == 'H3')
plot((H3$SEASON), H3$age, ylab = 'age', xlab = 'Season', type = 'p', main = 'H3N2')
age.lm = summary(lm(age ~ SEASON, data = H3)); age.lm
lines(H3$SEASON, age.lm$coefficients[1]+age.lm$coefficients[2]*H3$SEASON)
cor.test(H3$SEASON, H3$age, method = 'spearman')

plot((H3$SEASON), H3$BIRTHYEAR, col = 'blue', ylab = 'Birth year', xlab = 'Season')
BY.lm = summary(lm(BIRTHYEAR ~ SEASON, data = H3)); BY.lm
lines(H3$SEASON, BY.lm$coefficients[1]+BY.lm$coefficients[2]*H3$SEASON)
cor.test(H3$SEASON, H3$BIRTHYEAR, method = 'spearman')


## Young children might be causing the negative slop on birth year. Try exclusing children < 10
# Plot age vs. season
par(mfrow = c(2, 1))
H1 = subset(all.data, FLU_SUBTYPE == 'H1' & age > 9)
plot((H1$SEASON), H1$age, ylab = 'age', xlab = 'Season', type = 'p', main = 'H1N1')
age.lm = summary(lm(age ~ SEASON, data = H1)); age.lm
lines(H1$SEASON, age.lm$coefficients[1]+age.lm$coefficients[2]*H1$SEASON)
cor.test(H1$SEASON, H1$age, method = 'spearman')

plot((H1$SEASON), H1$BIRTHYEAR, col = 'blue', ylab = 'Birth year', xlab = 'Season')
BY.lm = summary(lm(BIRTHYEAR ~ SEASON, data = H1)); BY.lm
lines(H1$SEASON, BY.lm$coefficients[1]+BY.lm$coefficients[2]*H1$SEASON)
cor.test(H1$SEASON, H1$BIRTHYEAR, method = 'spearman')


par(mfrow = c(2, 1))
H3 = subset(all.data, FLU_SUBTYPE == 'H3' & age > 9)
plot((H3$SEASON), H3$age, ylab = 'age', xlab = 'Season', type = 'p', main = 'H3N2')
age.lm = summary(lm(age ~ SEASON, data = H3)); age.lm
lines(H3$SEASON, age.lm$coefficients[1]+age.lm$coefficients[2]*H3$SEASON)
cor.test(H3$SEASON, H3$age, method = 'spearman')

plot((H3$SEASON), H3$BIRTHYEAR, col = 'blue', ylab = 'Birth year', xlab = 'Season')
BY.lm = summary(lm(BIRTHYEAR ~ SEASON, data = H3)); BY.lm
lines(H3$SEASON, BY.lm$coefficients[1]+BY.lm$coefficients[2]*H3$SEASON)
cor.test(H3$SEASON, H3$BIRTHYEAR, method = 'spearman')



# Look at boxplots
par(mfrow = c(2, 1))
H1 = subset(all.data, FLU_SUBTYPE == 'H1')
plot(as.factor(H1$SEASON), H1$age, ylab = 'age', xlab = 'Season', type = 'p', main = 'H1N1')
age.lm = summary(lm(age ~ SEASON, data = H1)); age.lm
lines(H1$SEASON, age.lm$coefficients[1]+age.lm$coefficients[2]*H1$SEASON)
cor.test(H1$SEASON, H1$age, method = 'spearman')

plot(as.factor(H1$SEASON), H1$BIRTHYEAR, ylab = 'Birth year', xlab = 'Season')
BY.lm = summary(lm(BIRTHYEAR ~ SEASON, data = H1)); BY.lm
lines(H1$SEASON, BY.lm$coefficients[1]+BY.lm$coefficients[2]*H1$SEASON)
cor.test(H1$SEASON, H1$BIRTHYEAR, method = 'spearman')


par(mfrow = c(2, 1))
H3 = subset(all.data, FLU_SUBTYPE == 'H3')
plot(as.factor(H3$SEASON), H3$age, ylab = 'age', xlab = 'Season', type = 'p', main = 'H3N2')
age.lm = summary(lm(age ~ SEASON, data = H3)); age.lm
lines(H3$SEASON, age.lm$coefficients[1]+age.lm$coefficients[2]*H3$SEASON)
cor.test(H3$SEASON, H3$age, method = 'spearman')

plot(as.factor(H3$SEASON), H3$BIRTHYEAR, ylab = 'Birth year', xlab = 'Season')
BY.lm = summary(lm(BIRTHYEAR ~ SEASON, data = H3)); BY.lm
lines(H3$SEASON, BY.lm$coefficients[1]+BY.lm$coefficients[2]*H3$SEASON)
cor.test(H3$SEASON, H3$BIRTHYEAR, method = 'spearman')

