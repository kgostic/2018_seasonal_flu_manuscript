## Analyze Arizona severity data
rm(list = ls())
setwd('~/Dropbox/R/2017_seasonal_flu/')

# Import data
dat = read.csv('AZ_2016_Hospitalization.csv')
head(dat)

# import imprinting patterns
load('Seasonal_group_weights.RData')
wm.H1 = wm.H1['2016', ]
wo.H1 = wo.H1['2016', ]
wm.H3 = wm.H3['2016', ]
wo.H3 = wo.H3['2016', ]

# Add an age column to master data frame
dat$age = 2016-dat$BYEAR

# Add imprinting columns to master data frame
dat$p.g1.imp = round(wm.H1[as.character(dat$age)], 2)
dat$p.g2.imp = round(wm.H3[as.character(dat$age)], 2)

# Add a column saying whether the patient was admitted for 1 or more days
dat$admitted = I(dat$day_admit) != '0: no admit'

# Convert age gropus to factors
dat$Age_gp = factor(dat$Age_gp)

# indicate whether or not individuals are protected
#     --> >80% prob -> protected
#     --> <20% prob -> unprotected
#     --> 20<x<80% prob -> unknown
dat$H1.protected = 'unknown'
dat$H1.protected[which(dat$p.g1.imp > .8)] = 'protected'
dat$H1.protected[which(dat$p.g1.imp < .2)] = 'unprotected'

dat$H3.protected = 'unknown'
dat$H3.protected[which(dat$p.g2.imp > .8)] = 'protected'
dat$H3.protected[which(dat$p.g2.imp < .2)] = 'unprotected'

dat$H1.protected = factor(dat$H1.protected)
dat$H3.protected = factor(dat$H3.protected)



# Birth year categories
#     --> 1918-1956
#     --> 1957-1967
#     --> 1968-1977
#     --> 1978-2016
dat$BY.cat = NA
dat$BY.cat[which(dat$BYEAR >= 1918 & dat$BYEAR < 1957)] = '1918-1956'
dat$BY.cat[which(dat$BYEAR >= 1957 & dat$BYEAR < 1968)] = '1957-1967'
dat$BY.cat[which(dat$BYEAR >= 1968 & dat$BYEAR < 1977)] = '1968-1977'
dat$BY.cat[which(dat$BYEAR >= 1977 & dat$BYEAR < 2016)] = '1977-2016'

dat$BY.cat = factor(dat$BY.cat)


########################################################################
#########  TEST EFFECT OF IMPRINTING STATUS ON PROB ADMISSION  #########
########################################################################
# Explanatory factors include:
#  1. Age group ()
#  2. Imprinting status
#  3. Subtype
#
# Outcome is: admitted

# Separate H1N1 cases from H3N2
H1.dat = subset(dat, Subtype == 'H1N1')
H3.dat = subset(dat, Subtype == 'H3N2')

########## H1N1 analysis ##########
set.seed(31)

fit.0 = glm(I(admitted)~ns(age, df = 1)+H1.protected, data = H1.dat, family = binomial) # No knots = linear model
fit.1 = glm(I(admitted)~ns(age, df = 2)+H1.protected, data = H1.dat, family = binomial) # No knots = linear model
fit.2 = glm(I(admitted)~ns(age, df = 3)+H1.protected, data = H1.dat, family = binomial) # No knots = linear model
fit.3 = glm(I(admitted)~ns(age, df = 4)+H1.protected, data = H1.dat, family = binomial) # No knots = linear model
fit.4 = glm(I(admitted)~ns(age, df = 5)+H1.protected, data = H1.dat, family = binomial) # No knots = linear model
fit.5 = glm(I(admitted)~ns(age, df = 6)+H1.protected, data = H1.dat, family = binomial) # No knots = linear model


# Compare AIC
AICs = c(linear = fit.0$aic, k1 = fit.1$aic, k2 = fit.2$aic, k3 = fit.3$aic, k4 = fit.4$aic, k5 = fit.5$aic)
sort(AICs - min(AICs))
# One knot at 0.5 is best.



# 2. treat age as a continuous variable
library(gam)
par(mfrow = c(2, 3))
gam.fit = glm(I(admitted)~ns(age, df = 4)+H1.protected, data = H1.dat, family = binomial)
plot.gam(gam.fit, se = TRUE, col = 'blue')
title('H1N1')
plot(0:99, 1-wm.H1, col = 'blue', main = 'Prob unprotected\n by imprinting', type = 'l', lty = 2, xlab = 'age', ylab = 'prob unprotected')

gam.fit = glm(I(admitted)~H3.protected+ns(age, df = 4), data = H3.dat, family = binomial)
plot.gam(gam.fit, se = TRUE, col = 'red')
title('H3N2')
plot(0:99, 1-wm.H3, col = 'red', main = 'Prob unprotected\n by imprinting', type = 'l', lty = 2, xlab = 'age', ylab = 'prob unprotected')




## ICU admission
par(mfrow = c(2, 3))
gam.fit = glm(I(ICU)~H1.protected+ns(age, df = 3), data = H1.dat, family = binomial)
plot.gam(gam.fit, se = TRUE, col = 'blue')
title('H1N1')
plot(0:99, 1-wm.H1, col = 'blue', main = 'Prob unprotected\n by imprinting', type = 'l', lty = 2, xlab = 'age', ylab = 'prob unprotected')

gam.fit = glm(I(ICU)~H3.protected+ns(age, df = 3), data = H3.dat, family = binomial)
plot.gam(gam.fit, se = TRUE, col = 'red')
title('H3N2')
plot(0:99, 1-wm.H3, col = 'red', main = 'Prob unprotected\n by imprinting', type = 'l', lty = 2, xlab = 'age', ylab = 'prob unprotected')




# 
# ## Full model
# library(gam)
# gam.fit = glm(I(admitted)~ns(p.g1.imp, knots=c(1/3, 2/3))+Age_gp, data = H1.dat, family = binomial)
# par(mfrow = c(2, 2))
# plot.gam(gam.fit, se = TRUE, col = 'blue')
# title('H1N1')
# 
# gam.fit = glm(I(admitted)~ns(p.g2.imp, knots=c(1/3, 2/3))+Age_gp, data = H3.dat, family = binomial)
# plot.gam(gam.fit, se = TRUE, col = 'red')
# 
# 
# 
# ## Try another way:
table(H1.dat$H1.protected, H1.dat$admitted)

table(H1.dat$H1.protected, H1.dat$ICU)

table(H1.dat$Age_gp, H1.dat$admitted)



## BY Categories
## ICU admission
par(mfrow = c(2, 3))
gam.fit = glm(I(ICU)~BY.cat+ns(age, df = 3), data = H1.dat, family = binomial)
plot.gam(gam.fit, se = TRUE, col = 'blue')
title('H1N1')
plot(0:99, 1-wm.H1, col = 'blue', main = 'Prob unprotected\n by imprinting', type = 'l', lty = 2, xlab = 'age', ylab = 'prob unprotected')

gam.fit = glm(I(ICU)~BY.cat+ns(age, df = 3), data = H3.dat, family = binomial)
plot.gam(gam.fit, se = TRUE, col = 'red')
title('H3N2')
plot(0:99, 1-wm.H3, col = 'red', main = 'Prob unprotected\n by imprinting', type = 'l', lty = 2, xlab = 'age', ylab = 'prob unprotected')


