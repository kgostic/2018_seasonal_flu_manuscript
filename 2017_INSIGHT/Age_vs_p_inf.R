rm(list = ls())
setwd('~/Dropbox/R/2017_INSIGHT/')

source('Import_FLU002.R')

group.fun = function(N = 10462, grp.size = 10, min.ID = 0){
  n.grps = ceiling(N/grp.size)
  IDs = rep((1:n.grps)+min.ID, each = grp.size); IDs = IDs[1:N]
  IDs
}


ss = subset(dat.002, challenge == 'H3N2')
ss = ss[order(ss$age),]
IDs = group.fun(N = nrow(ss), grp.size = 20, min.ID = 0)
means = ages = vector('numeric', max(IDs))
for(ii in 1:max(IDs)){
  means[ii] = mean(ss$infected[IDs == ii])
  ages[ii] = mean(ss$age[IDs==ii])
}


## check your uninfected vs. infected labels


par(mfrow = c(3,1))

col2rgb('red')
mycol <- rgb(225, 165, 0, max = 255, alpha = 175)
c2 = rgb(225, 50, 0, max = 225, alpha = 150)

plot(ages, (means), cex = 1, col = mycol, pch = 16)

# get estimate and binomial CI for each age group
aa = sort(unique(ss$age))
CIlow = CIhigh = estimate = vector('numeric', length(aa))
for(ii in 1:length(aa)){
  valid = ss$infected[ss$age == aa[ii]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
}


tab = table(ss$age, ss$infected)
tab = tab[,c(2,1)]
xx = barplot(t(tab/rowSums(tab)), col = c(c2, mycol), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'Age vs.  freq. confirmed H3N2')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')


## Pool older age gropus
old = which(aa>=75)
means[old] = mean(means[old])
valid = ss$infected[ss$age >= 75]
bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
CIlow[old] = bb$conf.int[1]
CIhigh[old] = bb$conf.int[2]
estimate[old] = bb$estimate
tt = t(tab/rowSums(tab))
tt[1,old] = estimate[old[1]]
tt[2,old] = 1-estimate[old[1]]
xx = barplot(tt, col = c(c2, mycol), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'Age vs.  freq. confirmed H3N2')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 = CIlow, y1 = CIhigh, col = 'white')


# pool into 5 year age groups
summary(ss$age)
grps = matrix(17:96, ncol = 5, byrow = T)
bars = matrix(NA, nrow = 2, ncol = nrow(grps))
for(ii in 1:nrow(grps)){
  valid = ss[ss$age %in% grps[ii,],]
  bars[1,ii] = sum(valid$infected)
  bars[2,ii] = length(valid$infected) - sum(valid$infected)
}

estimate = apply(bars, 2, function(xx) binom.test(xx[1], n = sum(xx), conf.level = .95)$est)
CIlow = apply(bars, 2, function(xx) binom.test(xx[1], n = sum(xx), conf.level = .95)$conf.int[1])
CIhigh = apply(bars, 2, function(xx) binom.test(xx[1], n = sum(xx), conf.level = .95)$conf.int[2])


bars = t(t(bars)/colSums(bars))

xx = barplot(bars[c(1,2),], col = c(c2, mycol), border = NA, space = 0, xlab = 'age group', ylab = 'fraction', main = 'Age group vs.  freq. confirmed H3N2')
text(8, .1, 'Infected', col = 'white')
text(8, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 = CIlow, y1 = CIhigh, col = 'white')


rm(ss)
rm(tab)
rm(tt)


## Repeat for H1N1
ss = subset(dat.002, challenge == 'H1N1')
ss = ss[order(ss$age),]
IDs = group.fun(N = nrow(ss), grp.size = 20, min.ID = 0)
means = ages = vector('numeric', max(IDs))
for(ii in 1:max(IDs)){
  means[ii] = mean(ss$infected[IDs == ii])
  ages[ii] = mean(ss$age[IDs==ii])
}

col2rgb('dodgerblue')
mycol <- rgb(0, 150, 0, max = 255, alpha = 175)
c2 = rgb(10, 90, 150, max = 225, alpha = 150)

plot(ages, (means), cex = 1, col = mycol, pch = 16)

# get estimate and binomial CI for each age group
aa = sort(unique(ss$age))
CIlow = CIhigh = estimate = vector('numeric', length(aa))
for(ii in 1:length(aa)){
  valid = ss$infected[ss$age == aa[ii]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
}


tab = table(ss$age, ss$infected)
tab = tab[,c(2,1)]
xx = barplot(t(tab/rowSums(tab)), col = c(c2, mycol), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'Age vs.  freq. confirmed H1N1')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 = CIlow, y1 = CIhigh, col = 'white')



## Pool older age gropus
old = which(aa>=75)
means[old] = mean(means[old])
valid = ss$infected[ss$age >= 75]
bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
CIlow[old] = bb$conf.int[1]
CIhigh[old] = bb$conf.int[2]
estimate[old] = bb$estimate
tt = t(tab/rowSums(tab))
tt[1,old] = estimate[old[1]]
tt[2,old] = 1-estimate[old[1]]
xx = barplot(tt, col = c(c2, mycol), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'Age vs.  freq. confirmed H1N1')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 = CIlow, y1 = CIhigh, col = 'white')


# pool into 5 year age groups
summary(ss$age)
grps = matrix(17:96, ncol = 5, byrow = T)
bars = matrix(NA, nrow = 2, ncol = nrow(grps))
for(ii in 1:nrow(grps)){
  valid = ss[ss$age %in% grps[ii,],]
  bars[1,ii] = sum(valid$infected)
  bars[2,ii] = length(valid$infected) - sum(valid$infected)
}

estimate = apply(bars, 2, function(xx) binom.test(xx[1], n = sum(xx), conf.level = .95)$est)
CIlow = apply(bars, 2, function(xx) binom.test(xx[1], n = sum(xx), conf.level = .95)$conf.int[1])
CIhigh = apply(bars, 2, function(xx) binom.test(xx[1], n = sum(xx), conf.level = .95)$conf.int[2])


bars = t(t(bars)/colSums(bars))

xx = barplot(bars, col = c(c2, mycol), border = NA, space = 0, xlab = 'age group', ylab = 'fraction', main = 'Age group vs.  freq. confirmed H3N2')
text(8, .1, 'Infected', col = 'white')
text(8, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 = CIlow, y1 = CIhigh, col = 'white')












## Plot H1N1 and H3N2 infected counts across age groups
H1 = subset(dat.002, challenge == 'H1N1' & infected == 1)
H3 = subset(dat.002, challenge == 'H3N2' & infected == 1)

par(mfrow = c(2,1))
plot(as.numeric(names(table(H1$age))), table(H1$age), ylab = 'count', xlab = 'age', col = 'dodgerblue')
points(as.numeric(names(table(H3$age))), table(H3$age), col = 'firebrick1')

## now try overall age dist
H1c = subset(dat.002, challenge == 'H1N1')
H3c = subset(dat.002, challenge == 'H3N2')
plot(as.numeric(names(table(H1c$age))), table(H1c$age), ylab = 'count', xlab = 'age', col = 'dodgerblue')
points(as.numeric(names(table(H3c$age))), table(H3c$age), col = 'firebrick1')


barplot(table(H1c$infected, H1c$age), main = 'H1N1', ylab = 'counts')
tt = table(H1c$infected, H1c$age)
norm = (t(t(tt)/colSums(tt)))
barplot(norm, ylab = 'freq', xlab = 'age')



barplot(table(H3c$infected, H3c$age), main = 'H3N2', ylab = 'counts')
tt = table(H3c$infected, H3c$age)
norm = (t(t(tt)/colSums(tt)))
barplot(norm, ylab = 'freq', xlab = 'age')




