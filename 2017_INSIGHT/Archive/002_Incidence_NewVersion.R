## FLU 002 Incidence
## Prelim analysis

### transform the outcome variable using the logit link and make plots to explore the distribution of Yi with respect to each predictor.
setwd('~/Dropbox/R/2017_INSIGHT/')
rm(list = ls())
library(reshape)
source('Import_FLU002.R')
set.seed(43)
cols = c("#4477AA", "#DDCC77", "#CC6677")
cols = colors()[c(393, 464, 29)]
cols = colors()[c(506, 499, 475)]
# Pull out variables used in model
full.dat = dat.002
dat.002 = subset(dat.002, select = c('age', 'infected', 'anyvac', 'anydx', 'anyav', 'season', 'country', 'challenge', 'protected_grp', 'protected_sub', 'protected_N'))
which(apply(dat.002, MARGIN = 1,  function(xx) (any(is.na(xx)))))
# Exclude rows with any NAs
dat.002 = na.omit(dat.002)
# Exclude those with unknown vaccination status.
invalid = which(dat.002$anyvac == 2)
dat.002 = dat.002[-invalid, ]
# Exclude two observations with age < 18
invalid = which(dat.002$age == 17)
dat.002 = dat.002[-invalid,]

# Group countries with 30 or fewer valid cases into "other"
cx = names(which(table(dat.002$country) <= 30))
dat.002$country[dat.002$country %in% cx] = 'Other'

# Convert blocking variables to factors
dat.002$anyav = factor(dat.002$anyav)
dat.002$anyvac = factor(dat.002$anyvac)

dat.002$anydx = factor(dat.002$anydx)
dat.002$country = factor(dat.002$country)
dat.002$season = factor(dat.002$season)
dat.002$challenge = factor(dat.002$challenge)



par(mfrow = c(2,2))
####################
##  Vaccination
####################
vactab = table(dat.002$anyvac, dat.002$infected); colnames(vactab) = c('Inf=0', 'Inf=1')
xx = barplot(t(vactab/rowSums(vactab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'vaccinated', ylab = 'freq. infection (red)', main = 'a. 002 Fraction Infected')
CI0 = binom.test(x = vactab[1,2], n = rowSums(vactab)[1])$conf.int[1:2]
CI1 = binom.test(x = vactab[2,2], n = rowSums(vactab)[2])$conf.int[1:2]
segments(x0 = xx[1], y0 = CI0[1], y1 = CI0[2])
segments(x0 = xx[2], y0 = CI1[1], y1 = CI1[2])


####################
##  Antivirals
####################
avtab = table(dat.002$anyav, dat.002$infected); colnames(avtab) = c('Inf=0', 'Inf=1')
barplot(t(avtab/rowSums(avtab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'antivirals', ylab = 'freq. infection (red)', main = "b. 002 Fraction Infected")
CI0 = binom.test(x = avtab[1,2], n = rowSums(avtab)[1])$conf.int[1:2]
CI1 = binom.test(x = avtab[2,2], n = rowSums(avtab)[2])$conf.int[1:2]
segments(x0 = xx[1], y0 = CI0[1], y1 = CI0[2])
segments(x0 = xx[2], y0 = CI1[1], y1 = CI1[2])

####################
##  Underlying symptoms
####################
dxtab = table(dat.002$anydx, dat.002$infected); colnames(dxtab) = c('Inf=0', 'Inf=1')
barplot(t(dxtab/rowSums(dxtab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'underlying symptoms', ylab = 'freq. infection (red)', main = "c. 002 Fraction Infected")
CI0 = binom.test(x = dxtab[1,2], n = rowSums(dxtab)[1])$conf.int[1:2]
CI1 = binom.test(x = dxtab[2,2], n = rowSums(dxtab)[2])$conf.int[1:2]
segments(x0 = xx[1], y0 = CI0[1], y1 = CI0[2])
segments(x0 = xx[2], y0 = CI1[1], y1 = CI1[2])

####################
##  Challenge
####################
ctab = table(dat.002$challenge, dat.002$infected); colnames(ctab) = c('Inf=0', 'Inf=1')
barplot(t(ctab/rowSums(ctab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'Challenge', ylab = 'freq. infection (red)', main = "d. 002 Fraction Infected")
CI0 = binom.test(x = ctab[1,2], n = rowSums(ctab)[1])$conf.int[1:2]
CI1 = binom.test(x = ctab[2,2], n = rowSums(ctab)[2])$conf.int[1:2]
segments(x0 = xx[1], y0 = CI0[1], y1 = CI0[2])
segments(x0 = xx[2], y0 = CI1[1], y1 = CI1[2])


####################
##  Age
####################
# get estimate and binomial(link = "cloglog") CI for each age group
aa = sort(unique(dat.002$age))
CIlow = CIhigh = estimate = vector('numeric', length(aa))
for(ii in 1:length(aa)){
  valid = dat.002$infected[dat.002$age == aa[ii]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
}


tab = table(dat.002$age, dat.002$infected)
tab = tab[,c(2,1)]
xx = barplot(t(tab/rowSums(tab)), col = c('indianred3', 'orange'), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'a. Age vs. freq infection')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')
# Vaciance generally increases with age. Consider a var. stabilizing tns.
# #Exclude the two 17 year olds

pp = t(tab/rowSums(tab))[1,]
#plot(rownames(tab), pp)
lo = log(pp/(1-pp))
plot(rownames(tab), -lo, main = 'b. Age vs. neg log odds infection', xlab = 'age', ylab = 'negative log odds infection', col = 'firebrick1', cex = (rowSums(tab)-mean(rowSums(tab)))/max(rowSums(tab))*2+1)




####################
##  Age by challenge
####################
# get estimate and binomial(link = "cloglog") CI for each age group
H1 = dat.002[dat.002$challenge == 'H1N1', ]
aa = sort(unique(H1$age))
CIlow = CIhigh = estimate = vector('numeric', length(aa))
for(ii in 1:length(aa)){
  valid = H1$infected[H1$age == aa[ii]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
}

tab = table(H1$age, H1$infected)
tab = tab[,c(2,1)]
xx = barplot(t(tab/rowSums(tab)), col = c('slateblue', 'lightblue'), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'a. Age vs. freq infection\nH1N1 seasons')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')
# Vaciance generally increases with age. Consider a var. stabilizing tns.
# #Exclude the two 17 year olds
th1 = tab


H3 = dat.002[dat.002$challenge == 'H3N2', ]
aa = sort(unique(H3$age))
CIlow = CIhigh = estimate = vector('numeric', length(aa))
for(ii in 1:length(aa)){
  valid = H3$infected[H3$age == aa[ii]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
}

tab = table(H3$age, H3$infected)
tab = tab[,c(2,1)]
xx = barplot(t(tab/rowSums(tab)), col = c('indianred3', 'orange'), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'a. Age vs. freq infection\nH3N1 seasons')
text(35, .1, 'Infected', col = 'white')
text(35, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')
# Vaciance generally increases with age. Consider a var. stabilizing tns.
# #Exclude the two 17 year olds
th2 = tab

plot(as.numeric(rownames(th1)), th1[,1]/sum(th1[,1]), col = 'blue')
points(as.numeric(rownames(th2)), th2[,1]/sum(th2[,1]), col = 'red')


####################
##  Protection status
####################
aa = sort(unique(dat.002$protected_grp))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = dat.002$infected[dat.002$protected_grp > bins[ii,1] & dat.002$protected_grp <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
cc = (rowSums(tab)-mean(rowSums(tab)))/max(rowSums(tab))*2+1
summary(rowSums(tab))
tab = tab/rowSums(tab)


xx = barplot(t(tab), col = c('gray50', 'lightcyan3'), border = NA, space = 0, xlab = 'prob protected against challenge', ylab = 'fraction', main = 'c. Protection against H1N1 or H3N2', names.arg = binnames)
text(5, .1, 'Infected', col = 'white')
text(5, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')


pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'd. Protection against H1N1 or H3N2', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'slateblue', cex = cc)



pdf('INSIGHT_AgeDists_FullData.pdf')
par(mfrow = c(1,1), bg = 'black', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white', pch = 16)
plot(as.numeric(rownames(th1)), th1[,1]/sum(th1[,1]), col = 'dodgerblue', main = 'INSIGHT Outpatient\nCases from H1N1 or H3N2 dominated seasons only', xlab = 'Case age', ylab = 'Fraction of cases')
points(as.numeric(rownames(th2)), th2[,1]/sum(th2[,1]), col = 'red')
legend(70, .035, legend = c('H1N1', 'H3N2'), col = c('dodgerblue', 'red'), pch = 16, bty = 'n')
text(80, .038, paste(sum(th1[,1]), " H1N1\n", sum(th2[,1]), " H3N2"))


aa = subset(H1, infected == 1)
t1 = table(aa$age, aa$season)
t1 = t1[, which(colSums(t1)>100)]
bb = subset(H3, infected == 1)
t2 = table(bb$age, bb$season)
t2 = t2[, which(colSums(t2)>100)]
par(mfrow = c(1,1), bg = 'black', col.axis = 'white', col.lab = 'white', col.main = 'white', col = 'white', pch = 16)
xx=as.numeric(rownames(t1))
plot(xx, t1[,1]/sum(t1[,1]), col = 'dodgerblue', main = 'INSIGHT Outpatient\nCases from H1N1 or H3N2 dominated seasons only', xlab = 'Case age', ylab = 'Fraction of cases')
points(xx, t1[,2]/sum(t1[,2]), col = 'dodgerblue', pch = 2)
points(xx, t1[,3]/sum(t1[,3]), col = 'dodgerblue', pch = 3)
points(xx, t1[,4]/sum(t1[,4]), col = 'dodgerblue', pch = 4)
points(xx, t1[,5]/sum(t1[,5]), col = 'dodgerblue', pch = 5)
points(xx, t1[,6]/sum(t1[,6]), col = 'dodgerblue', pch = 6)
points(xx, t1[,7]/sum(t1[,7]), col = 'dodgerblue', pch = 7)


xx = rownames(t2)
points(xx, t2[,1]/sum(t2[,1]), col = 'firebrick1', pch = 16)
points(xx, t2[,2]/sum(t2[,2]), col = 'firebrick1', pch = 2)
points(xx, t2[,3]/sum(t2[,3]), col = 'firebrick1', pch = 3)
points(xx, t2[,4]/sum(t2[,4]), col = 'firebrick1', pch = 4)
points(xx, t2[,5]/sum(t2[,5]), col = 'firebrick1', pch = 5)
points(xx, t2[,6]/sum(t2[,6]), col = 'firebrick1', pch = 6)
points(xx, t2[,7]/sum(t2[,7]), col = 'firebrick1', pch = 7)
text(80, .038, paste(sum(t1), " H1N1\n", sum(t2), " H3N2"))
dev.off()





####################
##  Protection status by challenge
####################
H1 = dat.002[dat.002$challenge == 'H1N1', ]
aa = sort(unique(H1$protected_grp))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = H1$infected[H1$protected_grp > bins[ii,1] & H1$protected_grp <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
tab = tab/rowSums(tab)


xx = barplot(t(tab), col = c('navy', 'lightcyan3'), border = NA, space = 0, xlab = 'prob protected', ylab = 'fraction', main = 'p(protection at group level) vs.  freq. H1N1 infection', names.arg = binnames)
text(5, .1, 'Infected', col = 'white')
text(5, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')


pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')







H3 = dat.002[dat.002$challenge == 'H3N2', ]
aa = sort(unique(H3$protected_grp))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = H3$infected[H3$protected_grp > bins[ii,1] & H3$protected_grp <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
tab = tab/rowSums(tab)


xx = barplot(t(tab), col = c('darkred', 'orange'), border = NA, space = 0, xlab = 'prob protected', ylab = 'fraction', main = 'p(protection at group level) vs.  freq. H3N2 infection', names.arg = binnames)
text(5, .1, 'Infected', col = 'white')
text(5, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')



pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')



####################
##  Protection status (subtype)
####################
aa = sort(unique(dat.002$protected_sub))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = dat.002$infected[dat.002$protected_sub > bins[ii,1] & dat.002$protected_sub <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
tab = tab/rowSums(tab)


xx = barplot(t(tab), col = c('gray50', 'lightcyan3'), border = NA, space = 0, xlab = 'prob protected', ylab = 'fraction', main = 'p(protection at group level) vs.  freq. infection, any subtype', names.arg = binnames)
text(5, .1, 'Infected', col = 'white')
text(5, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')


pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')





####################
##  Protection status by challenge
####################
H1 = dat.002[dat.002$challenge == 'H1N1', ]
aa = sort(unique(H1$protected_sub))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = H1$infected[H1$protected_sub > bins[ii,1] & H1$protected_sub <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
tab = tab/rowSums(tab)


xx = barplot(t(tab), col = c('navy', 'lightcyan3'), border = NA, space = 0, xlab = 'prob protected', ylab = 'fraction', main = 'p(protection at group level) vs.  freq. H1N1 infection', names.arg = binnames)
text(5, .1, 'Infected', col = 'white')
text(5, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')


pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')







H3 = dat.002[dat.002$challenge == 'H3N2', ]
aa = sort(unique(H3$protected_sub))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = H3$infected[H3$protected_sub > bins[ii,1] & H3$protected_sub <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
tab = tab/rowSums(tab)


xx = barplot(t(tab), col = c('darkred', 'orange'), border = NA, space = 0, xlab = 'prob protected', ylab = 'fraction', main = 'p(protection at group level) vs.  freq. H3N2 infection', names.arg = binnames)
text(5, .1, 'Infected', col = 'white')
text(5, .9, 'Uninfected', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')



pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')



## Conclusion:
## Consider a variance stabilizing transform for age if possible
## Consider a non-linear relationship for age
## Protection seems to have a quadratic relationship wtih logit, and an interaction with challenge





####################
##  Full model using glmer. No additive components possible, and I'm not sure about variance structures.
####################
# library(lme4)
# m0 = glmer(infected ~ poly(protected_grp, 2)*challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m0)
# 
# # Test against a linear effect for protection
# m1 = glmer(infected ~poly(protected_grp, 1)*challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m1)
# 
# ## Problem: addition of age makes this model nearly unidentifiable. Add age wtih no interaction between protected and challenge to see that it has no clear effect.
# m2 = glmer(infected ~ protected_grp^2*challenge + age + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m2)
# 
# m3 = glmer(infected ~poly(protected_sub, 2)*challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m3)
# 
# # Test against a linear effect for protection
# m4 = glmer(infected ~poly(protected_sub, 1)*challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m4)
# 
# ## Test N protection
# m5 = glmer(infected ~poly(protected_N, 2)*challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m5)
# 
# # Test against a linear effect for protection
# m6 = glmer(infected ~poly(protected_N, 1)*challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m6)
# 
# ## Test interaction between N and group protection
# m7 = glmer(infected ~ protected_grp*protected_N+challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m7)
# 
# ## Test interaction between N and group protection
# m8 = glmer(infected ~ protected_sub*protected_N+challenge + anydx + anyav + anyvac*challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "logit"), control = glmerControl(optimizer = 'Nelder_Mead'))
# summary(m8)
# 
# 
# 
# summary(m0)$AIC
# summary(m1)$AIC
# summary(m2)$AIC
# summary(m3)$AIC  ## Subtype level, quadratic nearly as good
# summary(m4)$AIC
# summary(m5)$AIC  ## N imprinting, quadratic does best
# summary(m6)$AIC
# summary(m7)$AIC
# summary(m8)$AIC
# 
# 
# ##########################################################################
# # ## Gamm4 package: 
# # ##         PROS: Does not use PQL. Implements random effects in a straightforward way.
# # ##         CONS: Incredibly slow. Also, I'm not sure how weights work.
# # ##         Ultimate conclusion: Use method below for now, but consider using this method for publication. Results of the two packages are very similar.
# # library(gamm4)
# # m0b = gamm4(infected ~ s(protected_grp)+s(protected_grp, by = challenge) + challenge + age + anydx + anyav + anyvac, random = ~ (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# # summary(m0b$mer)
# # summary(m0b$gam)
# # anova(m0b$gam)
# ##########################################################################
# 
# 
# ## TRIED TO IMPLEMENT AGE-SPECIFIC WEIGHTS, BUT THE AIC BASICALLY DOUBLED. WHAT'S GOING ON?
# # ages = unique(dat.002$age)
# # agecounts = vector('numeric', nrow(dat.002))
# # for(ii in 1:length(ages)){
# #   ind = which(dat.002$age == ages[ii])
# #   agecounts[ind] = length(ind) # STore the number of observations of age i
# # }
# # 
# # m0c = gamm4(infected ~ s(protected_grp) + s(protected_grp, by = challenge) + age + anydx + anyav + anyvac, random = ~ (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"), weights = agecounts)

dev.off()
plot(dat.002$age, dat.002$protected_grp, main = 'Analyzed outpatient data', xlab = 'age', ylab = 'prob. Imprinting protection')

ss = subset(x = dat.002, subset = protected_grp < .1)
hist(ss$age, col = 'gray', main = 'Histogram of ages of people with <.1 prob')



## gam provides a possible alternative to gamm4, with a backdoor to random effects
# ##         PROS: Does not use PQL. Pretty fast.
# ##         CONS: I don't know how weights work. Also, uses a backdoor to implement random effects and doesn't offer a straightforward variance estimate.
# ##         Ultimate conclusion: Using this for now for speed, but may use gamm4 for pub. Results of the two are similar.
library(mgcv)
m0 = gam(infected ~ s(protected_grp, by = challenge, bs = 'cr') + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m0)

# Test against a linear effect for protection
m1 = gam(infected ~ protected_grp*challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m1)

m3= gam(infected ~  s(protected_sub, by = challenge, bs = 'cr') + challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m3)

# Test against a linear effect for protection
m4= gam(infected ~protected_sub*challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m4)

## Test N protection
m5= gam(infected ~s(protected_N, by = challenge, bs = 'cr') + challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m5)

# Test against a linear effect for protection
m6= gam(infected ~protected_N*challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m6)

## Test interaction between N and group protection
m7= gam(infected ~protected_N*protected_grp+challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m7)

## Test interaction between N and group protection
m8= gam(infected ~protected_N*protected_sub+challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m8)

AIC(m0)
AIC(m1)
AIC(m3)
AIC(m4)
AIC(m5)
AIC(m6)
AIC(m7)
AIC(m8)


#############################################
#  Best model contains nonlinear effect of N imprinting
#  Subtype-specific is almost as good
#  Effect is always significant for H1 and not H3
#   --> suggests that H1N1 childhood imprinting protects against H1, but H3N2 imprinting doesn't protect much against H3
#   --> similar to patterns in glycosylation analysis where lots of older people have strong H1N1 protection
############################################
summary(m5)
par(mfrow = c(2,3))
plot.gam(m5, main = "Nonlinear, N Imprinting")
par(mfrow = c(2,2))
gam.check(m5, type = 'pearson')
par(mfrow = c(2,3))
plot.gam(m3, main = "Nonlinear, sub Imprinting")
par(mfrow = c(2,2))
gam.check(m3, type = 'pearson')

AICsum = AIC(m0, m1, m3, m4, m5, m6, m7, m8)
AICsum$del.AIC = AICsum$AIC - min(AICsum$AIC)
AICsum[order(AICsum$del.AIC), ]

##### Remaining issues: 
## 1. Should use a var stabilizing transform on age (weights), but can't
## 3. How do I validate heterogeneity?
## 3. Nonlinear relationship between imp protection status may arise because of variant-specific patterns.
