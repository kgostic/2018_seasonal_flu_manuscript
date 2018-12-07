## FLU 002 prolonged
## Prelim analysis

rm(list = ls())
setwd('~/Dropbox/R/2017_INSIGHT/')
source('Import_FLU002.R')



# Notes on duration of symptoms:
#   A review by Carrat et al. (2007) of challenge studies in healthy human volunteers suggests that influenza symptoms typically resolve fully within 8-9 days of inoculation. This led us to classify patients with the following characteristics as having long-lasting symptoms:
#   
#   * Symptoms lasting 10 or more days prior to study enrollment (this implies that symptoms on day 10 were still serious enough to send patients to the doctor, and most likely lasted at least a few days longer.)
# * Symptoms not resolved on day 14 of enrollment (implies that symptoms lasted longer than 14 days, as most patients had symptoms for at least one day prior to enrollment.)

##############################
#  Define a binary variable: prolonged symptoms, yes or no?
#############################
dat.002$prolonged = as.numeric(dat.002$symdur >= 10 | dat.002$resolved == "0")
full.dat = dat.002 # Archive all columns

#####################################################
## Pull out variables used in model
## Drop cases in which the patient was NOT infected!
#####################################################
dat.002 = subset(dat.002, select = c('age', 'prolonged', 'anyvac', 'anydx', 'anyav', 'season', 'country', 'challenge', 'protected_grp', 'protected_sub', 'protected_N'), subset = infected == 1)
which(apply(dat.002, MARGIN = 1,  function(xx) (any(is.na(xx)))))
# Exclude rows with any NAs
dat.002 = na.omit(dat.002)
# Exclude those with unknown vaccination status.
invalid = which(dat.002$anyvac == 2)
dat.002 = dat.002[-invalid, ]
# Exclude two observations with very low age
invalid = which(dat.002$age == 17)
dat.002 = dat.002[-invalid,]




# Group countries with 30 or fewer valid cases into "other"
table(dat.002$country)
#dat.002$country[which(dat.002$country %in% c('Austria', 'Chile', 'Japan', 'UK', 'Portugal'))] = 'Other'

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
vactab = table(dat.002$anyvac, dat.002$prolonged); colnames(vactab) = c('Inf=0', 'Inf=1')
xx = barplot(t(vactab/rowSums(vactab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'vaccinated', ylab = 'freq. infection (red)', main = 'a. 002 Fraction w. Prolonged Symp')
# CI0 = binom.test(x = vactab[1,2], n = rowSums(vactab)[1])$conf.int[1:2]
# CI1 = binom.test(x = vactab[2,2], n = rowSums(vactab)[2])$conf.int[1:2]
# segments(x0 = xx[1], y0 = CI0[1], y1 = CI0[2])


####################
##  Antivirals
####################
avtab = table(dat.002$anyav, dat.002$prolonged); colnames(avtab) = c('Inf=0', 'Inf=1')
barplot(t(avtab/rowSums(avtab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'antivirals', ylab = 'freq. infection (red)', main = "b. 002 Fraction w. Prolonged Symp")

####################
##  Underlying symptoms
####################
dxtab = table(dat.002$anydx, dat.002$prolonged); colnames(dxtab) = c('Inf=0', 'Inf=1')
barplot(t(dxtab/rowSums(dxtab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'underlying symptoms', ylab = 'freq. infection (red)', main = "c. 002 Fraction w. Prolonged Symp.")

####################
##  Challenge
####################
ctab = table(dat.002$challenge, dat.002$prolonged); colnames(ctab) = c('Inf=0', 'Inf=1')
barplot(t(ctab/rowSums(ctab))[c(2,1),], border = NA, space = .1, col = c('indianred2', 'lightcyan3'), xlab = 'Challenge', ylab = 'freq. infection (red)', main = "d. 002 Fraction w. Prolonged Symp")


####################
##  Age
####################
# get estimate and binomial CI for each age group
aa = sort(unique(dat.002$age))
CIlow = CIhigh = estimate = vector('numeric', length(aa))
for(ii in 1:length(aa)){
  valid = dat.002$prolonged[dat.002$age == aa[ii]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
}


tab = table(dat.002$age, dat.002$prolonged)
tab = tab[,c(2,1)]
xx = barplot(t(tab/rowSums(tab)), col = c('indianred3', 'orange'), border = NA, space = 0, xlab = 'age', ylab = 'fraction', main = 'Age vs.  freq. prolonged, any subtype')
text(35, .1, 'prolonged', col = 'white')
text(35, .9, 'not prolonged', col = 'white')
points(xx, estimate, pch = '-', col = 'white')
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')
# Vaciance generally increases with age. Consider a var. stabilizing tns.
# #Exclude the two 17 year olds

pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot((as.numeric(rownames(tab))), -lo, main = 'Age', xlab = 'age', ylab = '-log odds prolonged', col = 'firebrick1', cex = (rowSums(tab)-mean(rowSums(tab)))/max(rowSums(tab))+1)


####################
##  Protection status
####################
aa = sort(unique(dat.002$protected_grp))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = dat.002$prolonged[dat.002$protected_grp > bins[ii,1] & dat.002$protected_grp <= bins[ii,2]]
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
text(5, .1, 'Prolonged', col = 'white')
text(5, .9, 'Not Prolonged', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')


pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), -lo, main = 'Protection status', xlab = 'p(protected)', ylab = '- log odds prolonged', col = 'lightcyan3', cex = (rowSums(tab)-mean(rowSums(tab)))/max(rowSums(tab))+1)



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
  valid = H1$prolonged[H1$protected_grp > bins[ii,1] & H1$protected_grp <= bins[ii,2]]
  bb = binom.test(x = sum(valid), n = length(valid), conf.level = .95)
  CIlow[ii] = bb$conf.int[1]
  CIhigh[ii] = bb$conf.int[2]
  estimate[ii] = bb$estimate
  inf[ii] = sum(valid)
  notinf[ii] = length(valid)-sum(valid)
}

tab = cbind(inf, notinf)
tab = tab/rowSums(tab)


pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')


xx = barplot(t(tab), col = c('navy', 'lightcyan3'), border = NA, space = 0, xlab = 'prob protected', ylab = 'fraction', main = 'p(protection at group level) vs.  freq. H1N1 infection', names.arg = binnames)
text(5, .1, 'Prolonged', col = 'white')
text(5, .9, 'Not prolonged', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')







H3 = dat.002[dat.002$challenge == 'H3N2', ]
aa = sort(unique(H3$protected_grp))
inf = notinf = CIlow = CIhigh = estimate = vector('numeric', 10)
bins = matrix(c(seq(0, .9, .1), seq(.1, 1, .1)), nrow = 10, ncol = 2)
binnames = apply(bins, 1, function(xx) paste(xx[1], xx[2], sep = '-'))
bins[1,1] = -.01
for(ii in 1:nrow(bins)){
  valid = H3$prolonged[H3$protected_grp > bins[ii,1] & H3$protected_grp <= bins[ii,2]]
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
text(5, .1, 'Prolonged', col = 'white')
text(5, .9, 'Not prolonged', col = 'white')
points(xx, estimate, pch = '-', col = 'white', cex = 1)
segments(x0 = xx, y0 =CIlow, y1 = CIhigh, col = 'white')

pp = t(tab/rowSums(tab))[1,]
lo = log(pp/(1-pp))
plot(rowMeans(bins), lo, main = 'Protection status', xlab = 'p(protected)', ylab = '(log odds)_age', col = 'lightcyan3')





## Conclusion:
## Expect to see a linear effect of age
## Expect to see no visible effect of protection status




# 
# ####################
# ##  Full model wtih no var transform for age vs. wtih var transform for age
# ####################
# library(lme4)
# library(nlme)
# m0 = glmer(prolonged ~ protected_grp + age + anyvac + anyav + anydx + challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m0)
# drop1(m0)
# 
# m00 = glmer(prolonged ~ protected_sub + age + anyvac + anyav + anydx + challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m00)
# 
# 
# m000 = glmer(prolonged ~ protected_N + age + anyvac + anyav + anydx + challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m000)
# 
# 
# # Drop protection
# m1 = glmer(prolonged ~ age + anyvac + anyav + anydx + challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m1)
# drop1(m1)
# 
# # Drop anyvac
# m2 = glmer(prolonged ~ age  + anyav + anydx + challenge + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m2)
# drop1(m2)
# 
# # Drop challenge
# m3 = glmer(prolonged ~ age  + anyav + anydx + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m3)
# drop1(m3)
# 
# # Drop anyav
# m4 = glmer(prolonged ~ age + anydx + (1|country)+(1|season), data = dat.002, family = binomial(link = "cloglog"))
# summary(m4)
# drop1(m4)
# 
# summary(m0)$AIC
# summary(m00)$AIC
# summary(m000)$AIC
# summary(m1)$AIC
# summary(m2)$AIC
# summary(m3)$AIC
# summary(m4)$AIC  # Least complex model wins!


## gam provides a possible alternative to gamm4, with a backdoor to random effects
# ##         PROS: Does not use PQL. Pretty fast.
# ##         CONS: I don't know how weights work. Also, uses a backdoor to implement random effects and doesn't offer a straightforward variance estimate.
# ##         Ultimate conclusion: Using this for now for speed, but may use gamm4 for pub. Results of the two are similar.
library(mgcv)
m0 = gam(prolonged ~ s(protected_grp, by = challenge, bs = 'cr') + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m0)
AIC(m0)

# Test against a linear effect for protection
m1 = gam(prolonged ~ protected_grp*challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m1)

m3= gam(prolonged ~  s(protected_sub, by = challenge, bs = 'cr') + challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m3)

# Test against a linear effect for protection
m4= gam(prolonged ~protected_sub*challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m4)

## Test N protection
m5= gam(prolonged ~s(protected_N, by = challenge, bs = 'cr') + challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m5)

# Test against a linear effect for protection
m6= gam(prolonged ~protected_N*challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m6)

## Test interaction between N and group protection
m7= gam(prolonged ~protected_N*protected_grp+challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m7)

## Test interaction between N and group protection
m8= gam(prolonged ~protected_N*protected_sub+challenge + age + anydx + anyav + anyvac*challenge + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m8)

AIC(m0)
AIC(m1)
AIC(m3)
AIC(m4)
AIC(m5)
AIC(m6)
AIC(m7)
AIC(m8)

## drop terms
## Test N protection
m5a = gam(prolonged ~s(protected_N, by = challenge, bs = 'cr') + age + anydx + anyav + anyvac + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m5a)
AIC(m5, m5a)

m5b = gam(prolonged ~s(protected_N, by = challenge, bs = 'cr') + age + anydx + s(country, bs = 're') + s(season, bs ="re"), data = dat.002, family = binomial(link = "logit"), select = TRUE)
summary(m5b)
AIC(m5, m5a, m5b)


summary(m5b)
par(mfrow = c(2,3))
plot.gam(m5b, main = "Nonlinear, N Imprinting")
par(mfrow = c(2,2))
gam.check(m5b, type = 'pearson')


##### Remaining issugs: 
## 1. Should use a var stabilizing transform on age (weights), but can't
## 3. How do I validate variance?
