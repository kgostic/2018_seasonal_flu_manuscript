library(ISLR)
Wage
head(Wage)


# Split the data in two and fit two separate models
w1 = Wage[1:1500, ]
w2 = Wage[1501:3000, ]
aa = glm(w1$wage ~ w1$age)
# w3 = w2; w3[] = NA
# w11 = rbind(w1, w3)
# aaa = lm(wage ~ age, data = w11)
bb = glm(w2$wage ~ w2$age)
# 
# ## Check NA grid
# colnames(w1) = paste(colnames(w1), rep(1, ncol(w1)), sep = '')
# colnames(w2) = paste(colnames(w2), rep(2, ncol(w2)), sep = '')
# w1empty = w1; w1empty[] = NA
# w2empty = w2; w2empty[] = NA
# www = rbind(cbind(w1, w2empty), cbind(w1empty, w2))
# www$wage = rowSums(cbind(www$wage1, www$wage2), na.rm = TRUE)
# lm(wage ~ age1 + age2, data = www)

# Add an indicator function to tell where the split happens
Wage$valid.1 = rep(c(1, 0), each = 1500)
Wage$valid.2 = rep(c(0, 1), each = 1500)
# Fit a combined lm and check that the results are the same as above
full1 = glm(Wage$wage ~ Wage$age + Wage$valid.2*Wage$age)
full2 = glm(Wage$wage ~ Wage$age*Wage$valid.1 + Wage$age)



## Summarize fits
caa = summary(aa)$coefficients
cbb = summary(bb)$coefficients
cf1 = summary(full1)$coefficients
cf2 = summary(full2)$coefficients

## Check that the coefficients and intercepts are the same in the individual and combined models.
aa.0 = caa[1, 1] # Null intercept, dataset 1
bb.0 = caa[2, 1] # Null coefficent, dataset 1
aa.f1 = cf1[1,1] # Full intercept for model where dataset 1 is baseline
bb.f1 = cf1[2, 1] # Full coefficient for model where dataset 1 is baseline
# check that the estimates are the same in the null model and full model 1
abs(aa.0 - aa.f1) < 1e-5 # Should be true
abs(bb.0 - bb.f1) < 1e-5
# check that the estimates are the same in the null model and the full model 2
aa.f2 = cf2[1,1] + cf2[3,1]
bb.f2 = cf2[2,1] + cf2[4,1]
# Add the base estimate wtih the interaction terms to get the estimates for dataset 2
abs(aa.0 - aa.f2) < 1e-5 # Should be true
abs(bb.0 - bb.f2) < 1e-5
# Convert standard errors
## VAR(A+B) = VAR(A) + VAR(B) + 2COV(A,B)
vcv1 = vcov(summary(full1))
vaa = vcv1[1,1] + vcv1[3,3] + 2*vcv1[1,3]; sqrt(vaa) # SE of the estimate for intercept, dataset 2
cf2[1, 2]
