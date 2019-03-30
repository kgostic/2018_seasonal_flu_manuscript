## Compare age dists by antigenic advance

## Clear memory
rm(list = ls())
#setwd('2017_AZ/')
library(lubridate)
library(ggplot2)
library(ggpubr)
library(reshape2)
setwd('~/Dropbox/R/2018_seasonal_flu/2017_AZ/')


## OUTPUTS


#######################################
## Load data, model inputs, and likelihood function
######################################
source('00-Inputs_multinomial.R')
## Load nextstrain antigenic advace data
H3N2_aa = read.table(file = 'nextstrain_data/nextstrain_staging_flu_seasonal_h3n2_ha_21y_metadata.tsv', sep = "\t", header = TRUE)
H1N1_aa_post2009 = read.table(file = 'nextstrain_data/nextstrain_flu_seasonal_h1n1pdm_ha_12y_metadata.tsv', sep = "\t", header = TRUE)
H1N1_aa_pre2009 = read.table(file = "nextstrain_data/elife-01914-fig3-data1-v1.tsv", sep = '\t', header = TRUE) 
H3N2_aa_bedford_elife = subset(H1N1_aa_pre2009, lineage == 'H3N2')
H1N1_aa_pre2009 = subset(H1N1_aa_pre2009, lineage == 'H1N1' & year > 2001) # Extract seasonal H1N1 cases from the relevant time period
H1N1_aa_pre2009$Full.date = as.character(H1N1_aa_pre2009$Full.date)
H1N1_aa_pre2009$Full.date = as.Date(H1N1_aa_pre2009$Full.date, format = '%m/%d/%y')
H1N1_aa_pre2009$decimal.date = decimal_date(H1N1_aa_pre2009$Full.date)
## Fill in entries with no clear date info using the raw year
replace = which(is.na(H1N1_aa_pre2009$decimal.date))
H1N1_aa_pre2009$decimal.date[replace] = H1N1_aa_pre2009$year[replace]


#####################################################################################################################
### Rescale estimates from bedford et al., eLife to match the magnitude of estimates from nextstrain (uses neher et al. methods)
##  H1N1 estimates are not available from the same time periods in both data sets, but H3N2 estimates from 1997-2011 are available in both data sets
##  Use paired H3N2 estimates to re-scale the Bedford et al. estimates to match nextstrain estimates
## First, find the mean antigenic location along dimension 1 (Bedford et al) per year
bedford_H3N2_yearly = sapply(1997:2011, function(xx){valid = H3N2_aa_bedford_elife$year == xx; mean(H3N2_aa_bedford_elife$ag1[valid])}); names(bedford_H3N2_yearly) = 1997:2011
##  Repeat for Neher et al
neher_H3N2_yearly = sapply(1997:2011, function(xx){valid = floor(H3N2_aa$Num.Date) == xx; mean(H3N2_aa$CTiter[valid], na.rm = TRUE)})

## Visualise the raw mean locations
par(mfrow = c(1,2))
plot(1997:2011, bedford_H3N2_yearly-min(bedford_H3N2_yearly), xlab = 'year', ylab = 'antigenic location', main = 'Raw estimates')
points(1997:2011, neher_H3N2_yearly-min(neher_H3N2_yearly), col = 'red')
legend('topleft', c('bedford', 'neher'), col = c('black', 'red'), pch = 1)
## Find the scaling factor that standardizes estimates to span the same range
scale_factor = diff(range(neher_H3N2_yearly))/diff(range(bedford_H3N2_yearly))
## Rescale and visualize the rescaled points
plot(1997:2011, (bedford_H3N2_yearly-min(bedford_H3N2_yearly))*scale_factor, xlab = 'year', ylab = 'antigenic location', main = 'Rescaled Bedford estimates; Neher unchanged')
points(1997:2011, neher_H3N2_yearly-min(neher_H3N2_yearly), col = 'red')

## Rescale the bedford H1N1 estimates
H1N1_aa_pre2009$ag1 = H1N1_aa_pre2009$ag1*scale_factor



#####################################################################################################################
## CALCULATE ANTIGENIC ADVANCE FOR FOUR CATEGORIES:
##      1. H3N2 (all years)
##      2. H1N1, pre-2009 (re-scaled to match nextstrain/neher estimates)
##      3. 2009 pandemic H1N1 (assume antigenic advance = 0 since this is a new strain)
##      4. H1N1 post-2009 pandemic
####################################################################################################################

############################################# H3N2 ############################################# 
## - Separate isolates into influenza seasons
## - Find the mean antigenic location per season
## - Find the season-to-season differnce between means to measure average antigenic advance per year
## Establish a baseline using data from nextstrain
lowest = floor(min(H3N2_aa$Num.Date, na.rm = TRUE))
## Our definition is that the NH influenza sesaon begins in week 40
##  There are 52.14 weeks in a year
##  NH flu season begins in decimal week 40/52.14 = 0.77
##  Separate all specimens into NH flu seasons, starting with the 1997-1998 season
##  Find average divergence from strains that circulated in 1997, prior to week 40
##  Then find the difference between season-to-season divergence
yrs = 1997:max(floor(H3N2_aa$Num.Date), na.rm = TRUE) # Years of interest
CTiter.raw = CTiterSub.raw = numeric(length(yrs)) # Initialize raw divergence
names(CTiter.raw) = names(CTiterSub.raw) = paste(yrs, yrs+1, sep = "-")
baseline = colMeans(subset(H3N2_aa, Num.Date < 1997.77, select = c('CTiter', 'CTiterSub')))
for(yy in 1:length(yrs)){
  valid = which(H3N2_aa$Num.Date >= yrs[yy]+.77 & H3N2_aa$Num.Date < yrs[yy]+1.77) # Extract all the sample indices from a given NH season (week 40-week 39)
  CTiter.raw[yy] = mean(H3N2_aa[valid, 'CTiter'])
  CTiterSub.raw[yy] = mean(H3N2_aa[valid, 'CTiterSub'])
}
plot(yrs, CTiter.raw)
points(yrs, CTiterSub.raw)
## Find the season-to-season difference
CTiter.H3N2 = c('1997-1998' = 0, diff(CTiter.raw))
CTiterSub.H3N2 = c('1997-1998' = 0, diff(CTiterSub.raw))
plot(seq(1997.5, 2018.5, by = 1), CTiter.H3N2, col = 'red', main = 'H3N2 CTiter'); abline(h = 0)
plot(seq(1997.5, 2018.5, by = 1), CTiterSub.H3N2, col = 'red', main = 'H3N2 CTiterSub'); abline(h = 0)


############################################# Post-pandemic H1N1 ############################################# 
## ---------------------  aa measured on same scale as H3N2 using nextstrain data ----------------------------
## Establish a baseline
yrs = 2009:2018 # Years of interest
CTiter.raw = CTiterSub.raw = numeric(length(yrs)) # Initialize raw divergence
names(CTiter.raw) = names(CTiterSub.raw) = paste(yrs, yrs+1, sep = "-")
baseline = colMeans(subset(H1N1_aa_post2009, Num.Date < lowest+.77, select = c('CTiter', 'CTiterSub')))
for(yy in 1:length(yrs)){
  valid = which(H1N1_aa_post2009$Num.Date >= yrs[yy]+.77 & H1N1_aa_post2009$Num.Date < yrs[yy]+1.77) # Extract all the sample indices from a given NH season (week 40-week 39)
  CTiter.raw[yy] = mean(H1N1_aa_post2009[valid, 'CTiter'])
  CTiterSub.raw[yy] = mean(H1N1_aa_post2009[valid, 'CTiterSub'])
}
plot(yrs, CTiter.raw, col = 'blue', main = 'H1N1 post-pandemic CTiter')
points(yrs, CTiterSub.raw, col = 'red', main = 'H1N1 post-pandemic CTiterSub')
## Find the season-to-season difference
CTiter.H1N1 = c('2009-2010' = 0, diff(CTiter.raw))
CTiterSub.H1N1 = c('2009-2010' = 0, diff(CTiterSub.raw))
plot(seq(2009.5, 2018.5, by = 1), CTiter.H1N1, col = 'blue'); abline(h = 0)
points(seq(2009.5, 2018.5, by = 1), CTiterSub.H1N1, col = 'red'); abline(h = 0)



############################################# pandemic H1N1 ############################################# 
## -------------------------- Only one year of circulation, aa = NA ------------------------------


############################################# pre-pandemic H1N1 ############################################# 
## ----------------- aa originally measured on a different scale than post-pandemic H1N1 and H3N2-----------------------
##  antigenic location esimates were rescled above to match post-pandemic and H3N2 estimates

yrs = 2002:2008 # Years of interest
aassn.raw = numeric(length(yrs)) # Initialize raw divergence
names(aassn.raw) =  paste(yrs, yrs+1, sep = "-")
baseline = colMeans(subset(H1N1_aa_pre2009, decimal.date < lowest+.77, select = c('ag1')))
for(yy in 1:length(yrs)){
  valid = which(H1N1_aa_pre2009$decimal.date >= yrs[yy]+.77 & H1N1_aa_pre2009$decimal.date < yrs[yy]+1.77) # Extract all the sample indices from a given NH season (week 40-week 39)
  aassn.raw[yy] = mean(H1N1_aa_pre2009[valid, 'ag1'])
}

plot(yrs, aassn.raw, col = 'blue')

## Find the season-to-season difference
aassn = c('2002-2003' = 0, diff(aassn.raw))
plot(seq(2002.5, 2008.5, by = 1), aassn, col = 'blue'); abline(h = 0)












#######################################
## Reformat H3N2 data for plotting
######################################
H3.master = H3.master[-(1:2),] # Drop seassons before 200203 season, the first year with corresponding data on antigenic advance.
## Remove data from earlier seasons
num.season = as.numeric(gsub(rownames(H3.master), pattern = '(\\d{4})\\d{2}', replacement = '\\1'))+1 ## extract the second year of each sesaon. Use this to convert birth years to ages
age.mat = t(sapply(num.season, FUN = function(xx){xx-2015:1918})) # Entires correspond to the age of each entry in H3.master
## Enter case counts for ages 0-85 into each season starting with 02-03
H3.dat = matrix(NA, nrow = nrow(H3.master), ncol = 86, dimnames = list(rownames(H3.master), 0:85))
for(ii in 1:nrow(H3.dat)){
  H3.dat[ii, ] = H3.master[ii, which(age.mat[ii,] %in% 0:85)]
}
## Get rid of seasons with fewer than 50 observations
H3.dat = H3.dat[rowSums(H3.dat)>=70, ]
## Calculte frequencies for plotting
H3.freq = H3.dat/rowSums(H3.dat)
H3.freq_cumulative = t(apply(H3.freq, 1, cumsum))

## Melt the age-specific frequencies into a data frame, with antigenic advance ('CTiter') as ID vars
H3.summary =   melt(as.data.frame(cbind(freq = H3.freq,
                                        'CTiter' = CTiter.H3N2[gsub(rownames(H3.dat), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")], 
                                        'CTiterSub' = CTiterSub.H3N2[gsub(rownames(H3.dat), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")],
                                        'season' = rownames(H3.freq))), id.vars = c('CTiter', 'CTiterSub', 'season'), variable.name = 'age', value.name = 'frequency')

## Add cumulative frequency
H3.summary$c.freq =  melt(as.data.frame(cbind(c.freq = H3.freq_cumulative,
                                              'season' = rownames(H3.dat))), 
                          id.vars = 'season', variable_name = 'age', value_name = 'c.freq')$value
## Add counts                       
H3.summary$count =  melt(as.data.frame(cbind(count = H3.dat,
                                             'season' = rownames(H3.dat))), 
                         id.vars = 'season', variable_name = 'age', value_name = 'count')$value


## Define a function to reformat data to create histograms
age_bins = function(age_tab){
  out = cbind(rowSums(age_tab[,as.character(0:10)]),
        rowSums(age_tab[,as.character(11:25)]),
        rowSums(age_tab[,as.character(26:40)]),
        rowSums(age_tab[,as.character(41:55)]),
        rowSums(age_tab[,as.character(56:70)]),
        rowSums(age_tab[,as.character(71:85)]))
  out = out/rowSums(out)
  colnames(out) = c('0-10', '11-25', '26-40', '41-55', '56-70', '71-75')
  as.data.frame(out)
}
## Bin case counts into broad age groups
H3.age.bins = cbind(age_bins(H3.dat), 'CTiter' = CTiter.H3N2[gsub(rownames(H3.dat), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")], 'season' = rownames(H3.dat))










#######################################
## Refromat post-pandemic H1N1 for plotting
## Because we have antigenic advance data on the same scale as H3N2 data, we can plot these on the same panel
######################################
H1.master = H1.master[-(1:2),] # Drop seassons before 200203 season, the first year with corresponding data on antigenic advance.
## Drop seasons in which 0 cases were observed
H1.master = H1.master[-c(2,3), ]
## Extract post-pandemic cases
H1.post.pandemic = H1.master[5:9, ]


## Remove data from earlier seasons
num.season = as.numeric(gsub(rownames(H1.post.pandemic), pattern = '(\\d{4})\\d{2}', replacement = '\\1'))+1 ## extract the second year of each sesaon. Use this to convert birth years to ages
age.mat = t(sapply(num.season, FUN = function(xx){xx-2015:1918})) # Entires correspond to the age of each entry in H3.master
## Enter case counts for ages 0-85 into each season starting with 02-03
H1.dat = matrix(NA, nrow = nrow(H1.post.pandemic), ncol = 86, dimnames = list(rownames(H1.post.pandemic), 0:85))
for(ii in 1:nrow(H1.dat)){
  H1.dat[ii, ] = H1.post.pandemic[ii, which(age.mat[ii,] %in% 0:85)]
}
## Get rid of seasons with fewer than 50 observations
H1.dat = H1.dat[rowSums(H1.dat)>=70, ]
## Calculte frequencies for plotting
H1.freq = H1.dat/rowSums(H1.dat)
H1.freq_cumulative = t(apply(H1.freq, 1, cumsum))

## Melt the age-specific frequencies into a data frame, with antigenic advance ('CTiter') as ID vars
H1.summary =   melt(as.data.frame(cbind(freq = H1.freq,
                                        'CTiter' = CTiter.H1N1[gsub(rownames(H1.dat), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")], 
                                        'CTiterSub' = CTiterSub.H1N1[gsub(rownames(H1.dat), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")],
                                        'season' = rownames(H1.freq))), id.vars = c('CTiter', 'CTiterSub', 'season'), variable.name = 'age', value.name = 'frequency')

## Add cumulative frequency
H1.summary$c.freq =  melt(as.data.frame(cbind(c.freq = H1.freq_cumulative,
                                              'season' = rownames(H1.dat))), 
                          id.vars = 'season', variable_name = 'age', value_name = 'c.freq')$value
## Add counts                       
H1.summary$count =  melt(as.data.frame(cbind(count = H1.dat,
                                             'season' = rownames(H1.dat))), 
                         id.vars = 'season', variable_name = 'age', value_name = 'count')$value


H1.age.bins = cbind(age_bins(H1.dat), 'CTiter' = CTiter.H1N1[gsub(rownames(H1.dat), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")], 'season' = rownames(H1.dat))





#######################################
## Refromat pre-pandemic H1N1 for plotting
## Because we have antigenic advance data on the same scale as H3N2 data, we can plot these on the same panel
######################################
## Extract post-pandemic cases
H1.pre.pandemic = H1.master[1:4, ]


## Remove data from earlier seasons
num.season = as.numeric(gsub(rownames(H1.pre.pandemic), pattern = '(\\d{4})\\d{2}', replacement = '\\1'))+1 ## extract the second year of each sesaon. Use this to convert birth years to ages
age.mat = t(sapply(num.season, FUN = function(xx){xx-2015:1918})) # Entires correspond to the age of each entry in H3.master
## Enter case counts for ages 0-85 into each season starting with 02-03
H1.dat.pre = matrix(NA, nrow = nrow(H1.pre.pandemic), ncol = 86, dimnames = list(rownames(H1.pre.pandemic), 0:85))
for(ii in 1:nrow(H1.dat.pre)){
  H1.dat.pre[ii, ] = H1.pre.pandemic[ii, which(age.mat[ii,] %in% 0:85)]
}
## Get rid of seasons with fewer than 50 observations
H1.dat.pre = H1.dat.pre[rowSums(H1.dat.pre)>=70, ]
## Calculte frequencies for plotting
H1.freq.pre = H1.dat.pre/rowSums(H1.dat.pre)
H1.freq.pre_cumulative = t(apply(H1.freq.pre, 1, cumsum))

## Melt the age-specific frequencies into a data frame, with antigenic advance ('CTiter') as ID vars
H1.summary.pre =   melt(as.data.frame(cbind(freq = H1.freq.pre,
                     'CTiter' = aassn[gsub(rownames(H1.dat.pre), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")], 
                     'CTiterSub' = aassn[gsub(rownames(H1.dat.pre), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")],
                     'season' = rownames(H1.freq.pre))), id.vars = c('CTiter', 'CTiterSub', 'season'), variable.name = 'age', value.name = 'frequency')

## Add cumulative frequency
H1.summary.pre$c.freq =  melt(as.data.frame(cbind(c.freq = H1.freq.pre_cumulative,
                                    'season' = rownames(H1.dat.pre))), 
                                    id.vars = 'season', variable_name = 'age', value_name = 'c.freq')$value

## Add counts                       
H1.summary.pre$count =  melt(as.data.frame(cbind(count = H1.dat.pre,
                                              'season' = rownames(H1.dat.pre))), 
                          id.vars = 'season', variable_name = 'age', value_name = 'count')$value


H1.pre.age.bins = cbind(age_bins(H1.dat.pre), 'CTiter' = aassn[gsub(rownames(H1.dat.pre), pattern = "(\\d{4})(\\d{2})", replacement = "\\1-20\\2")], season = rownames(H1.dat.pre))




H1.age.bins = melt(H1.age.bins, id.vars = c('CTiter', 'season')); H1.age.bins$type = 'H1N1_post_2009'
H3.age.bins = melt(H3.age.bins, id.vars = c('CTiter', 'season')); H3.age.bins$type = 'H3N2'
H1.pre.age.bins = melt(H1.pre.age.bins, id.vars = c('CTiter', 'season')); H1.pre.age.bins$type = 'H1N1_seasonal'
full.age.bins = rbind(H1.age.bins, H3.age.bins, H1.pre.age.bins); full.age.bins$seasontype = paste(full.age.bins$season, full.age.bins$type, sep = '_')

## Reorder seasons according to CTiter
useasontype = unique(full.age.bins$seasontype)
uaa = unique(full.age.bins$CTiter)
full.age.bins$seasontype = factor(full.age.bins$seasontype,
                                  levels = useasontype[order(uaa)])

ggplot()+
  geom_bar(stat = 'identity', data = full.age.bins, aes(x = variable, y = value, fill = CTiter, group = seasontype, color = type), position = 'dodge')+
  scale_fill_viridis_c(option = 'plasma', na.value = 'gray')+
  scale_discrete_manual(values = c('black', 'gray', 'white'), aesthetics = 'color')

## Calculate pearson correlation coefficients for each variable
## Only calculate for H3N2 because there are too few data points for other types
get.cor = function(xx){
  valid = subset(H3.age.bins, variable == xx)
  out = cor.test(valid$CTiter, valid$value, method = 'pearson')
  c(r = as.numeric(out$estimate), p = as.numeric(out$p.value), variable = xx)
}
cor.df = as.data.frame(t(sapply(unique(H3.age.bins$variable), FUN = get.cor)))
cor.df$variable = unique(H3.age.bins$variable)
cor.df$label = paste('r=', round(cor.df$r,2), '  p=', round(cor.df$p,2), sep = '')
cor.df$type = 'H3N2'
cor.df$CTiter = NA

  ggplot()+
  facet_wrap(~variable) +
  geom_smooth(data = full.age.bins, aes(x = CTiter, y = value, group = type, color = type), method = 'lm', na.rm = TRUE, lwd = .5, lty = 2, se = FALSE) +
  geom_point(data = full.age.bins, aes(x = CTiter, y = value, fill = CTiter, color = type, shape = type)) +
  scale_color_manual(values = c('gray50', 'gray50', 'black')) +
  scale_fill_viridis_c(option = 'plasma', na.value = 'gray') +
  scale_shape_manual(values = c(21,22,23)) +
  theme_bw() +
  geom_label(data = cor.df, aes(x = 0, y = .6, label = label), hjust = 0) +
  xlab('Antigenic advance, relative to previous season') +
  ylab('Fraction of cases in age group')
  
  
anova_like_plot =  ggplot()+
    facet_wrap(~variable) +
    #geom_smooth(data = full.age.bins, aes(x = CTiter, y = value, group = type, color = type), method = 'lm', na.rm = TRUE, lwd = .5, lty = 3, se = FALSE) +
    geom_smooth(data = subset(full.age.bins, type == 'H3N2' & variable == '0-10'), aes(x = CTiter, y = value, group = type, color = type), method = 'lm', na.rm = TRUE, lwd = 1, lty = 2, se = FALSE, show.legend = FALSE) +
    geom_smooth(data = subset(full.age.bins, type == 'H3N2'), aes(x = CTiter, y = value, group = type, color = type), method = 'lm', na.rm = TRUE, lwd = .5, lty = 2, se = FALSE, show.legend = FALSE) +
    geom_point(data = full.age.bins, aes(x = CTiter, y = value, color = type, shape = type)) +
    theme_bw() +
    geom_label(data = cor.df, aes(x = 0, y = .65, label = label), hjust = 0, color = 4) +
    xlab('Antigenic advance, relative to previous season') +
    ylab('Fraction of cases in age group')

ggsave(filename = '../figures/Antigenic_advance_scatter.pdf', height = 7, width = 9)




remove_levels = function(df.in){
  df.in$age = as.numeric(as.character(df.in$age))
  df.in$freqquency = as.numeric(as.character(df.in$frequency))
  df.in$c.freq = as.numeric(as.character(df.in$c.freq))
  df.in$count = as.numeric(as.character(df.in$count))
  df.in$CTiter = as.numeric(as.character(df.in$CTiter))
  df.in$CTiterSub = as.numeric(as.character(df.in$CTiterSub))
  df.in
}



## Combine H3N2 data and post-pandemic H1N1 data
H1.summary= remove_levels(H1.summary)
H1.summary.pre = remove_levels(H1.summary.pre)
H3.summary = remove_levels(H3.summary)
H1.summary$type = 'H1N1pandemic'
H1.summary.pre$type = 'H1N1seasonal'
H3.summary$type = 'H3N2'
fulldat = rbind(H1.summary, H3.summary, H1.summary.pre)
fulldat$age = as.numeric(as.character(fulldat$age))
fulldat$seasontype = paste(fulldat$type, fulldat$season, sep = '_')


ggplot(data = fulldat, aes(x = age, y = c.freq)) + 
  geom_line(aes(color = CTiter, group = seasontype, linetype = type)) +
  #geom_point(aes(color = CTiter, group = seasontype, shape = type)) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma",  na.value = 'gray') 


ggplot(data = H1.summary, aes(x = age, y = c.freq)) + 
  geom_line(aes(color = CTiter, group = season), linetype = 1) +
  #geom_point(aes(color = CTiter, group = seasontype, shape = type)) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") 

ggplot(data = H1.summary.pre, aes(x = age, y = c.freq)) + 
  geom_line(aes(color = CTiter, group = season), linetype = 3) +
  #geom_point(aes(color = CTiter, group = seasontype, shape = type)) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") 

ggplot(data = H3.summary, aes(x = age, y = c.freq)) + 
  geom_line(aes(color = CTiter, group = season), linetype = 2) +
  #geom_point(aes(color = CTiter, group = seasontype, shape = type)) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") 








### CALCULATE KS STATISTIC FOR EACH PAIR OF FLU SEASONS
## Write a functiont to convert age table into long form
long_dat = function(short_dat){
  rep(0:85, short_dat)
}
## Write a function to calculate the KS statistic given long data inputs
KS_stat = function(ld1, ld2){
  ks.test(x = ld1, y = ld2)$statistic
}
## Write a function to calculate the KS pvalue given long data inputs
KS_pval = function(ld1, ld2){
  ks.test(x = ld1, y = ld2)$p.value
}



## Make a list of all data from after the 2009 pandemic (same aa scores)
H1.list.post = apply(X = H1.dat, MARGIN = 1, FUN = long_dat)
H1.list.pre = apply(H1.dat.pre, MARGIN = 1, FUN = long_dat)
H3.list = apply(H3.dat, MARGIN = 1, FUN = long_dat)
dat.table = rbind(H1.freq[,1:86], H3.freq[,1:86])
long.dat.list = c(H1.list.post, H3.list)
dat.types = rep(c('H1pdm', 'H3'), c(length(H1.list.post), length(H3.list)))
aa = c(H1.freq$CTiter, H3.freq$CTiter); names(aa )= paste(dat.types, names(long.dat.list), sep = '_') # Get antigenic advance scores


## Calculate the KS stat for each possible combination of seasons
nn = length(long.dat.list)
KS.mat = SSE.mat = matrix(NA, nrow = nn, ncol = nn, dimnames = list(paste(dat.types, names(long.dat.list), sep = '_'), paste(dat.types, names(long.dat.list), sep = '_')))
for(cc in 1:nn){
  for(rr in 1:nn){
    KS.mat[rr, cc] = ks.test(long.dat.list[[rr]], long.dat.list[[cc]])$statistic
    SSE.mat[rr, cc] = sum((dat.table[rr, ] - dat.table[cc,])^2)
  }
}


KS = melt(KS.mat)
KS$aa1 = sapply(KS$X1, FUN = function(xx) aa[which(names(aa) == xx)])
KS$aa2 = sapply(KS$X2, FUN = function(xx) aa[which(names(aa) == xx)])
KS$del.aa = abs(KS$aa1-KS$aa2)


SSE = melt(SSE.mat)
SSE$aa1 = sapply(SSE$X1, FUN = function(xx) aa[which(names(aa) == xx)])
SSE$aa2 = sapply(SSE$X2, FUN = function(xx) aa[which(names(aa) == xx)])
SSE$del.aa = abs(SSE$aa1-SSE$aa2)




ggplot()+
  geom_raster(data = SSE, aes(x = X1, y = X2, fill = (value)), )+
  scale_fill_viridis_c(option = 'plasma', limits = c(1e-8, max(SSE$value)+.01), na.value = 'gray')




