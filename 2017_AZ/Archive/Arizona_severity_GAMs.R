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
# For each imprinting status, use natural splines
# Use AIC to determine how many knots are best

# folds = sample(1:nrow(H1.dat), size = nrow(H1.dat), replace = FALSE)
# begin.folds = seq(1, 1000, by = 104)
# end.folds = begin.folds+104; end.folds[10] = nrow(H1.dat)

# # Set test and training dat
# test.dat = H1.dat[folds[begin.folds[1]:end.folds[1]],]
# training.dat = H1.dat[-folds[begin.folds[1]:end.folds[1]],]
# Fit a natural spline with various knots
fit.0 = glm(I(admitted)~ns(p.g1.imp, knots = NULL), data = H1.dat, family = binomial) # No knots = linear model
fit.1 = glm(I(admitted)~ns(p.g1.imp, knots = c(.5)), data = H1.dat, family = binomial) # One knot at prob = .5
fit.2 = glm(I(admitted)~ns(p.g1.imp, knots = c(1/3, 2/3)), data = H1.dat, family = binomial) # Two knots at prob = 1/3, 2/3

# Compare AIC
AICs = c(linear = fit.0$aic, k1 = fit.1$aic, k2 = fit.2$aic)
sort(AICs - min(AICs))
# One knot at 0.5 is best.

# Predict logit
preds.0 = predict(fit.0, newdata = list(p.g1.imp = seq(0, 1, .01)))
preds.1 = predict(fit.1, newdata = data.frame(p.g1.imp = seq(0, 1, .01)))
preds.2 = predict(fit.2, newdata = data.frame(p.g1.imp = seq(0, 1, .01)))

# Convert to probs
p.fit.0 = exp(preds.0)/(1+exp(preds.0))
p.fit.1 = exp(preds.1)/(1+exp(preds.1))
p.fit.2 = exp(preds.2)/(1+exp(preds.2))

plot(seq(0, 1, .01), p.fit.0)

# Draw outcomes
outcomes.0 = sapply(p.fit.0, function(pp) rbinom(1, 1, pp))
outcomes.1 = sapply(p.fit.1, function(pp) rbinom(1, 1, pp))
outcomes.2 = sapply(p.fit.2, function(pp) rbinom(1, 1, pp))


