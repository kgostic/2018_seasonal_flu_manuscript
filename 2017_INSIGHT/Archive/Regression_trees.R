## Build boosted regression trees for INSIGHT data
setwd('~/Dropbox/R/2017_INSIGHT/')
rm(list = ls())
source('Import_FLU002.R')
set.seed(2)
head(dat.002)
dat.002$anyav = factor(dat.002$anyav)
dat.002$anyvac = factor(dat.002$anyvac)
dat.002$anydx = factor(dat.002$anydx)
dat.002$season = factor(dat.002$season)
dat.002$country = factor(dat.002$country)


pdf('regressionTrees.pdf')


#######################################################
## Predict the outcome H1N1 vs. no PCR confirmed flu
## ----------------------------------------------------
##   Extract relevant columns
##  Split into a training and a test set
H1N1.incidence = subset(dat.002, subset = !is.na(H1N1),  select = c('H1N1', 'age', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'birth.year'))
train = sample(1:nrow(H1N1.incidence), nrow(H1N1.incidence)*7/8) # Train on 7/8 of the data
test = H1N1.incidence[-train, ]

## Boosted regression tree
library(gbm)
boost.H1N1 = gbm(H1N1 ~., data = H1N1.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H1N1)
plot(boost.H1N1, i = 'birth.year', col = 'blue')
abline(v = c(1957, 1968, 1977), lty = 2, col = 'black')
text(c(1957, 1968, 1977)-2, y = -1.7, c('1957', '1968', '1977'), srt = 90, col = 'black')
segments(x0 = 1918, x1 = 1957, y0= -1.63, col = 'dodgerblue')
text(1930, -1.62, 'H1N1', col = 'dodgerblue')
segments(x0 = 1957, x1 = 1968, y0= -1.63, col = 'blue')
text(1963, -1.62, 'H2N2', col = 'blue')
segments(x0 = 1968, x1 = 1977, y0= -1.63, col = 'firebrick1')
text(1972, -1.62, 'H3N2', col = 'firebrick1')
segments(x0 = 1977, x1 = 2015, y0= -1.63, col = 'purple')
text(1990, -1.62, 'Mixed', col = 'purple')

predvars = colnames(H1N1.incidence)[-c(1, 8)]
for(ii in 1:6){
  plot(boost.H1N1, i = predvars[ii])
}








#######################################################
## Predict the outcome H3N2 vs. no PCR confirmed flu
## ----------------------------------------------------
##   Extract relevant columns
##  Split into a training and a test set
H3N2.incidence = subset(dat.002, subset = !is.na(H3N2),  select = c('H3N2', 'birth.year', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'age'))
train = sample(1:nrow(H3N2.incidence), nrow(H3N2.incidence)*7/8) # Train on 7/8 of the data
test = H3N2.incidence[-train, ]


## Boosted regression tree
boost.H3N2 = gbm(H3N2 ~., data = H3N2.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H3N2)
plot(boost.H3N2, i = 'birth.year', col = 'blue')
abline(v = c(1957, 1968, 1977), lty = 2, col = 'black')
text(c(1957, 1968, 1977)-2, y = -1.32, c('1957', '1968', '1977'), srt = 90, col = 'black')
segments(x0 = 1918, x1 = 1957, y0= -1.275, col = 'dodgerblue')
text(1930, -1.27, 'H1N1', col = 'dodgerblue')
segments(x0 = 1957, x1 = 1968, y0= -1.275, col = 'blue')
text(1963, -1.27, 'H2N2', col = 'blue')
segments(x0 = 1968, x1 = 1977, y0= -1.275, col = 'firebrick1')
text(1972, -1.27, 'H3N2', col = 'firebrick1')
segments(x0 = 1977, x1 = 2015, y0= -1.275, col = 'purple')
text(1990, -1.27, 'Mixed', col = 'purple')
predvars = colnames(H3N2.incidence)[-c(1, 2)]
for(ii in 1:6){
  plot(boost.H3N2, i = predvars[ii])
}





## Build boosted regression trees for INSIGHT data
setwd('~/Dropbox/R/2017_INSIGHT/')
rm(list = ls())
source('Import_FLU002.R')
set.seed(2)
head(dat.002)
dat.002$anyav = factor(dat.002$anyav)
dat.002$anyvac = factor(dat.002$anyvac)
dat.002$anydx = factor(dat.002$anydx)
dat.002$season = factor(dat.002$season)
dat.002$country = factor(dat.002$country)








#######################################################
## Predict the outcome H1N1 vs. no PCR confirmed flu, using pg1 protection
## ----------------------------------------------------
##   Extract relevant columns
##  Split into a training and a test set
H1N1.incidence = subset(dat.002, subset = !is.na(H1N1),  select = c('H1N1', 'age', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'p.g1.protection'))
train = sample(1:nrow(H1N1.incidence), nrow(H1N1.incidence)*7/8) # Train on 7/8 of the data
test = H1N1.incidence[-train, ]

## Boosted regression tree
library(gbm)
boost.H1N1 = gbm(H1N1 ~., data = H1N1.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H1N1)
predvars = colnames(H1N1.incidence)[-c(1)]
for(ii in 1:7){
  plot(boost.H1N1, i = predvars[ii])
}








#######################################################
## Predict the outcome H3N2 vs. no PCR confirmed flu
## ----------------------------------------------------
##   Extract relevant columns
##  Split into a training and a test set
H3N2.incidence = subset(dat.002, subset = !is.na(H3N2),  select = c('H3N2', 'p.g2.protection', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'age'))
train = sample(1:nrow(H3N2.incidence), nrow(H3N2.incidence)*7/8) # Train on 7/8 of the data
test = H3N2.incidence[-train, ]



## Boosted regression tree
boost.H3N2 = gbm(H3N2 ~., data = H3N2.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H3N2)
predvars = colnames(H3N2.incidence)[-c(1)]
for(ii in 1:7){
  plot(boost.H3N2, i = predvars[ii])
}














#######################################################
## Predict the outcome H1N1 vs. no PCR confirmed flu, using pg1 protection
## ----------------------------------------------------
##   Extract relevant columns
##  Split into a training and a test set
H1N1.incidence = subset(dat.002, subset = !is.na(H1N1),  select = c('H1N1', 'age', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'p.g1.protection'))
train = sample(1:nrow(H1N1.incidence), nrow(H1N1.incidence)*7/8) # Train on 7/8 of the data
test = H1N1.incidence[-train, ]

## Boosted regression tree
library(gbm)
boost.H1N1 = gbm(H1N1 ~., data = H1N1.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H1N1)
predvars = colnames(H1N1.incidence)[-c(1)]
for(ii in 1:7){
  plot(boost.H1N1, i = predvars[ii], return.grid = T)
}








#######################################################
## Predict the outcome ICU vs. no ICU 003
## ----------------------------------------------------
##   Extract relevant columns
##  Split into a training and a test set
source('Import_FLU003.R')
H1N1.incidence = subset(dat.003, subset = H1N1 == 1,  select = c('anyicu', 'p.g2.protection', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'age'))
H1N1.incidence$anyav = factor(H1N1.incidence$anyav)
H1N1.incidence$anyvac = factor(H1N1.incidence$anyvac)
H1N1.incidence$anydx = factor(H1N1.incidence$anydx)
H1N1.incidence$season = factor(H1N1.incidence$season)
H1N1.incidence$country = factor(H1N1.incidence$country)
train = sample(1:nrow(H1N1.incidence), nrow(H1N1.incidence)*7/8) # Train on 7/8 of the data
test = H1N1.incidence[-train, ]



## Boosted regression tree
boost.H1N1 = gbm(anyicu ~., data = H1N1.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H1N1)
predvars = colnames(H1N1.incidence)[-c(1)]
for(ii in 1:7){
  plot(boost.H1N1, i = predvars[ii])
}









H3N2.incidence = subset(dat.003, subset = H3N2 == 1,  select = c('anyicu', 'p.g2.protection', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'age'))
H3N2.incidence$anyav = factor(H3N2.incidence$anyav)
H3N2.incidence$anyvac = factor(H3N2.incidence$anyvac)
H3N2.incidence$anydx = factor(H3N2.incidence$anydx)
H3N2.incidence$season = factor(H3N2.incidence$season)
H3N2.incidence$country = factor(H3N2.incidence$country)
train = sample(1:nrow(H3N2.incidence), nrow(H3N2.incidence)*7/8) # Train on 7/8 of the data
test = H3N2.incidence[-train, ]



## Boosted regression tree
boost.H3N2 = gbm(anyicu ~., data = H3N2.incidence[train, ], distribution = 'bernoulli', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.H3N2)
predvars = colnames(H3N2.incidence)[-c(1)]
for(ii in 1:7){
  plot(boost.H3N2, i = predvars[ii])
}





###########  3-way classification H1N1, H3N2 or 0
clean.dat = subset(dat.003, subset = flutype %in% c(1, 2, 4),  select = c('flutype', 'p.g1.protection', 'anyvac', 'anyav', 'anydx', 'country', 'season', 'age'))
clean.dat$anyav = factor(clean.dat$anyav)
clean.dat$anyvac = factor(clean.dat$anyvac)
clean.dat$anydx = factor(clean.dat$anydx)
clean.dat$season = factor(clean.dat$season)
clean.dat$country = factor(clean.dat$country)
clean.dat$flutype = factor(clean.dat$flutype)
train = sample(1:nrow(clean.dat), nrow(clean.dat)*7/8) # Train on 7/8 of the data
test = H3N2.incidence[-train, ]
## Boosted regression tree
boost.flutype = gbm(flutype ~., data = clean.dat[train, ], distribution = 'multinomial', n.trees = 5000, interaction.depth = 4)


## Partial dependence plots
par(mfrow = c(3,3), mar = c(5, 5, 1, 1), las = 2, cex.lab = .7)
summary(boost.flutype)
predvars = boost.flutype$var.names
for(ii in 1:7){
  plot(boost.flutype, i = predvars[ii])
}


