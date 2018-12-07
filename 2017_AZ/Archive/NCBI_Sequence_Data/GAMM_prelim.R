#### Visualize two kinds of sequence data and test GAMMS
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')

## Load NCBI data
H1_NCBI.raw = read.csv('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/NCBI_H1_Data.csv', stringsAsFactors = FALSE)
H3_NCBI.raw = read.csv('~/Dropbox/R/2017_seasonal_flu/NCBI_Sequence_Data/NCBI_H3_Data.csv', stringsAsFactors = FALSE)

##
H1_NCBI = subset(H1_NCBI.raw, subset = Protein == 'HA')
H3_NCBI = subset(H3_NCBI.raw, subset = Protein == 'HA')

## Load Arizona data
H1_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H1.csv')
H3_AZ = read.csv('~/Dropbox/R/2017_seasonal_flu/AZ_H3.csv')


## Plot H1 and H3, color by year
base.cols = c('lightseagreen', 'deeppink4', 'olivedrab1', 'red', 'yellow', 'dodgerblue', 'magenta', 'turquoise', 'goldenrod', 'purple', 'chocolate1', 'limegreen', 'darkblue')
cols = 'start'
for(ii in 1:13){
  ramp <- colorRamp(c(base.cols[ii], "white"))
  cols = c(cols, rgb( ramp(seq(0, .5, length = 2)), max = 255))
}
cols = cols[-1]

#reord = rev(c(1, 5, 9, 13, 17, 21, c(1, 5, 9, 13, 17, 21)+1, c(1, 5, 9, 13, 17, 21)+2, c(1, 5, 9, 13, 17, 21)+3))
reord = (c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26))
cols = cols[reord]
names(cols) = 2015:1990

yrs = as.character(1990:2015)


format.table = function(raw.data){
  rows = sort(unique(raw.data$Year), decreasing = TRUE)
  columns = 2015:1918
  outs = matrix(NA, nrow = length(rows), ncol = length(columns), dimnames = list(rows, columns))
  for(rr in 1:length(rows)){
    sub = subset(raw.data, Year == rows[rr])
    for(cc in 1:length(columns)){
    outs[rr, cc] = sum(sub$BirthYear == columns[cc])
    }
  }
  outs
}
    
    

####### Birth year
par(mfrow = c(2,2), mar = c(3, 3, 3, 1), bg = 'black', fg = 'white', col.main = 'white', col.axis = 'white', col.lab = 'white')
## H1 AZ
dt = format.table(H1_AZ); use.cols = which(names(cols) %in% rownames(dt))
xx = barplot(dt, col = cols[use.cols], border = NA, main = 'H1N1 Arizona', space = 0, xaxt = 'n')
axis(side = 1, at = xx[seq(1, 100, by = 15)], labels = seq(2015, 1920, by = -15))
mtext('Birth Year', side = 1, cex = .8, line = 2)
mtext('Case count', side = 2, line = 2, cex = .8)
mtext(paste('N=', sum(dt)), side = 3, line = 0, col = 'white', cex = .8)

## H1 NCBI
dt = format.table(H1_NCBI); use.cols = which(names(cols) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H1N1 NCBI', space = 0, xaxt = 'n')
axis(side = 1, at = xx[seq(1, 100, by = 15)], labels = seq(2015, 1920, by = -15))
mtext('Birth Year', side = 1, cex = .8, line = 2)
mtext(paste('N=', sum(dt)), side = 3, line = 0, col = 'white', cex = .8)

## H3 AZ
dt = format.table(H3_AZ); use.cols = which(names(cols) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H3N2 Arizona', space = 0, xaxt = 'n')
axis(side = 1, at = xx[seq(1, 100, by = 15)], labels = seq(2015, 1920, by = -15))
mtext('Birth Year', side = 1, cex = .8, line = 2)
mtext('Case count', side = 2, line = 2, cex = .8)
mtext(paste('N=', sum(dt)), side = 3, line = 0, col = 'white', cex = .8)

## H3 NCBI
dt = format.table(H3_NCBI); use.cols = which(names(cols) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H3N2 NCBI', space = 0, xaxt = 'n')
axis(side = 1, at = xx[seq(1, 100, by = 15)], labels = seq(2015, 1920, by = -15))
mtext('Birth Year', side = 1, cex = .8, line = 2)
mtext(paste('N=', sum(dt)), side = 3, line = 0, col = 'white', cex = .8)

plot.new(); par(mfrow = c(1,1))
plot.window(xlim = c(0,1), ylim = c(0,1))
legend('left', legend = as.character(1990:2002), fill = cols[as.character(1990:2002)], bty = 'n')
legend('right', legend = as.character(2003:2014), fill = cols[as.character(2003:2014)], bty = 'n')


## Density
par(mar = c(4, 4, 2, 1))
dens_H1 = density(H1_AZ$BirthYear, from = 1918, to = 2015)
dens_H3 = density(H3_AZ$BirthYear, from = 1918, to = 2015)
plot(dens_H1$x, dens_H1$y, type = 'l', col = 'dodgerblue', lwd = 2, ylab = 'Density', xlab = 'Birth year')
lines(dens_H3$x, dens_H3$y, col = 'red')
legend('topleft', c('H1N1\n', 'H3N2\n'), lty = 1, col = c('dodgerblue', 'red'), bty = 'n')


par(mar = c(4, 4, 2, 1))
dens_H1 = density(H1_AZ$Age, from = 0, to = 100)
dens_H3 = density(H3_AZ$Age, from = 0, to = 100)
plot(dens_H1$x, dens_H1$y, type = 'l', col = 'blue', lwd = 2, ylab = 'Density', xlab = 'Age')
lines(dens_H3$x, dens_H3$y, col = 'red')
legend(70, .027, c('H1N1\n', 'H3N2\n'), lty = 1, col = c('blue', 'red'), bty = 'n')





### Age
## H1 AZ
par(mfrow = c(2,2))
dt = table(H1_AZ$Year, H1_AZ$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H1 Arizona')

## H1 NCBI
dt = table(H1_NCBI$Year, H1_NCBI$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H1 NCBI')

## H3 AZ
dt = table(H3_AZ$Year, H3_AZ$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H3 Arizona')

## H3 NCBI
dt = table(H3_NCBI$Year, H3_NCBI$Age); use.cols = which(as.character(1990:2014) %in% rownames(dt))
barplot(dt, col = cols[use.cols], border = NA, main = 'H3 NCBI')






#### Both age patterns look the same
# AZ is from a single state
# NCBI is worldwide
library(gam)
library(splines)




## Not sure if this is statistically valid, but make a table of all possible birth years and ages

## Start with AZ H1 data!!
counts = table(H1_AZ$year, H1_AZ$BIRTHYEAR)/rowSums(table(H1_AZ$year, H1_AZ$BIRTHYEAR))
bys = rep(as.numeric(colnames(counts)), each = nrow(counts))
ages = as.vector(sapply(as.numeric(colnames(counts)), function(by) as.numeric(rownames(counts)) - by))
ages[which(ages < 0)] = NA
counts = as.vector(counts)
GAM.H1 = data.frame(props = counts, birth.year = bys, age = ages)
GAM.H1 = na.omit(GAM.H1)




## Determine the optimal number of df for age
## test 1:20 degrees of freedom
NN = nrow(GAM.H1) ## Do 11-fold cross validation
nn = ceiling(NN/11)
rand.ind = matrix(sample(x = 1:NN, size = NN, replace = FALSE), ncol = nn, nrow = NN/nn)


age.df = function(df.in){
fits = vector('list', 11)
rmse = vector('numeric', 11)
for(ff in 1:11){
  test.inds = rand.ind[ff, ]
  training = GAM.H1[-test.inds, ]
  test = GAM.H1[test.inds, ]
  fits[[ff]] = gam(props ~ ns(age, df = df.in), data = training)
  preds = predict(fits[[ff]], newdata=test)
  rmse[ff] = sum( sqrt((test$props - preds)^2 ))
}
sum(rmse)
}

age.cross.validation = sapply(1:20, age.df); names(age.cross.validation) = paste(1:20, 'df', sep = '')
plot(1:20, age.cross.validation, xlab = 'df', ylab = 'root mean square test error', type = 'b')
## use 8 or 9 df for age!     



by.df = function(df.in){
  fits = vector('list', 11)
  rmse = vector('numeric', 11)
  for(ff in 1:11){
    test.inds = rand.ind[ff, ]
    training = GAM.H1[-test.inds, ]
    test = GAM.H1[test.inds, ]
    fits[[ff]] = gam(props ~ ns(birth.year, df = df.in), data = training)
    preds = predict(fits[[ff]], newdata=test)
    rmse[ff] = sum( sqrt((test$props - preds)^2 ))
  }
  sum(rmse)
}
by.cross.validation = sapply(1:20, by.df)
plot(1:20, by.cross.validation, xlab = 'df', ylab = 'root mean square test error', type = 'b')
## use 4 df for age! 

rm(gm.by, gm.a, gm.by.a)


## First fit natural splines with various df for AGE
gm.by.a = gam(props ~ ns(birth.year, df = 4)+ ns(age, df = 9), data = GAM.H1) # N.cases ~ birth year
par(mfrow=c(2,1))
plot(gm.by.a, se=TRUE,col="blue ")
summary(gm.by.a)

## age only
gm.a = gam(props ~ ns(age, df = 9), data = GAM.H1) # N.cases ~ birth year
par(mfrow=c(1,1))
plot(gm.a, se=TRUE,col="blue ")
summary(gm.a)

## by only
gm.by = gam(props ~ ns(birth.year, df = 4), data = GAM.H1) # N.cases ~ birth year
par(mfrow=c(1,1))
plot(gm.by, se=TRUE,col="blue ")
summary(gm.by)

## Test full model vs. age only
anova(gm.a, gm.by.a, test = 'F')









## Now with AZ H3 data
counts = table(H3_AZ$year, H3_AZ$BIRTHYEAR)/rowSums(table(H3_AZ$year, H3_AZ$BIRTHYEAR))
bys = rep(as.numeric(colnames(counts)), each = nrow(counts))
ages = as.vector(sapply(as.numeric(colnames(counts)), function(by) as.numeric(rownames(counts)) - by))
ages[which(ages < 0)] = NA
counts = as.vector(counts)
GAM.H3 = data.frame(props = counts, birth.year = bys, age = ages)
GAM.H3 = na.omit(GAM.H3)




## Determine the optimal number of df for age
## test 1:20 degrees of freedom
NN = nrow(GAM.H3) ## Do 11-fold cross validation
nn = round(NN/10)
rand.ind = matrix(sample(x = 1:NN, size = NN-1, replace = FALSE), ncol = nn, nrow = 10)


age.df = function(df.in){
  fits = vector('list', 10)
  rmse = vector('numeric', 10)
  for(ff in 1:10){
    test.inds = rand.ind[ff, ]
    training = GAM.H3[-test.inds, ]
    test = GAM.H3[test.inds, ]
    fits[[ff]] = gam(props ~ ns(age, df = df.in), data = training)
    preds = predict(fits[[ff]], newdata=test)
    rmse[ff] = sum( sqrt((test$props - preds)^2 ))
  }
  sum(rmse)
}

age.cross.validation = sapply(1:20, age.df); names(age.cross.validation) = paste(1:20, 'df', sep = '')
plot(1:20, age.cross.validation, xlab = 'df', ylab = 'root mean square test error', type = 'b')
## use 11 df for age!     



by.df = function(df.in){
  fits = vector('list', 10)
  rmse = vector('numeric', 10)
  for(ff in 1:10){
    test.inds = rand.ind[ff, ]
    training = GAM.H3[-test.inds, ]
    test = GAM.H3[test.inds, ]
    fits[[ff]] = gam(props ~ ns(birth.year, df = df.in), data = training)
    preds = predict(fits[[ff]], newdata=test)
    rmse[ff] = sum( sqrt((test$props - preds)^2 ))
  }
  sum(rmse)
}

by.cross.validation = sapply(1:20, by.df)
plot(1:20, by.cross.validation, xlab = 'df', ylab = 'root mean square test error', type = 'b')
## use 9 df for age! 


## First fit natural splines with various df for AGE
gm.by.a = gam(props ~ ns(birth.year, df = 9)+ ns(age, df = 11), data = GAM.H3) # N.cases ~ birth year
par(mfrow=c(2,1))
plot(gm.by.a, se=TRUE,col="blue ")
summary(gm.by.a)

## age only
gm.a = gam(props ~ ns(age, df = 11), data = GAM.H3) # N.cases ~ birth year
par(mfrow=c(1,1))
plot(gm.a, se=TRUE,col="blue ")
summary(gm.a)

## by only
gm.by = gam(props ~ ns(birth.year, df = 9), data = GAM.H3) # N.cases ~ birth year
par(mfrow=c(1,1))
plot(gm.by, se=TRUE,col="blue ")
summary(gm.by)

## Test full model vs. age only
anova(gm.a, gm.by.a, test = 'F')


