---
title: "003 death"
author: "Katie Gostic"
date: "11/1/2017"
output: html_document
---


Notes on data formatting from Deborah Wentworth:

* symdur is symptom duration (days) prior to enrollment
* date info is for enrollment date
* fluvac6 and fluvac12 are flu vaccination in last 6 or last 12 months
* anyvac combines values from fluvac6 and fluvac12, whichever was collected for a given patient is reported.
* anyav is use of antivirals at any time
* anydx is any diagnosis in the medical history section
* flutype is subtype: 1=H1N1, 2=H3N2, 3=flu B, 4=PCR negative, 5=Inf A, subtype unknown, and 6=coinfection with multiple strains
* resolved is symptom resolution within 14 days
* country code is the 3 letter code for country of enrollment - think the abbreviations will be obvious


#### Visualize the relationship between age and prob of H1N1 infection
![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)



#### Visualize the relationship between age and prob of H3N2 infection
![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

#### Visualize the relationship between country and prolonged symptoms

H1N1

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

H3N2

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)


#### Visualize the relationship between season and prob infection

H1N1

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

H3N2

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)



### Clean data for fitting

```
## 
## Argentina Australia   Austria   Belgium     Chile     China   Denmark 
##        51        59         5         3         1        16        70 
##   Germany    Greece    Norway      Peru    Poland Singapore     Spain 
##        38        82         7        17         6         6        43 
##  Thailand        UK       USA 
##        64       130       111
```

```
## [1] 709
```

```
## 
## Argentina Australia   Belgium     China   Denmark   Germany    Greece 
##        93       113         3         7        30        21        20 
##      Peru Singapore     Spain  Thailand        UK       USA 
##         9        21        28       121        91       168
```

```
## [1] 725
```


### First fit a model where the only independent variable is age
Treat vaccination, antivial use, underlying symptoms, country and season as blocking variables with fixed effects

```r
library(splines)
## H1N1 fit to country, season and medical history only (no effect of age or imprinting)
fit0 = glm(dth ~ anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit0)
```

```
## 
## Call:
## glm(formula = dth ~ anyav + anydx + country + H3.valid * season + 
##     H3.valid * anyvac, family = binomial, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.0524  -0.2850  -0.1900  -0.0869   3.1654  
## 
## Coefficients: (1 not defined because of singularities)
##                            Estimate Std. Error z value Pr(>|z|)   
## (Intercept)                -4.00169    1.56213  -2.562  0.01042 * 
## anyav1                      0.47994    0.44233   1.085  0.27791   
## anydx1                      3.03302    1.03106   2.942  0.00326 **
## countryAustralia           -2.65391    0.97819  -2.713  0.00667 **
## countryDenmark             -1.43832    1.11520  -1.290  0.19714   
## countryGermany             -1.09444    1.18677  -0.922  0.35643   
## countryGreece              -1.80297    1.14047  -1.581  0.11390   
## countryOther               -1.14632    0.99881  -1.148  0.25110   
## countrySpain               -2.91714    1.46892  -1.986  0.04704 * 
## countryThailand            -3.50636    1.15903  -3.025  0.00248 **
## countryUK                  -2.06963    1.06516  -1.943  0.05201 . 
## countryUSA                 -2.64430    1.10405  -2.395  0.01662 * 
## H3.valid1                  13.17769 2957.38283   0.004  0.99644   
## seasonNH.10.11              0.05913    0.66405   0.089  0.92904   
## seasonNH.11.12              2.83153    1.39436   2.031  0.04228 * 
## seasonNH.12.13            -15.11671 2171.04675  -0.007  0.99444   
## seasonNH.13.14             -0.97438    1.16058  -0.840  0.40115   
## seasonNH.14.15            -15.22509 1609.33673  -0.009  0.99245   
## seasonNH.15.16              0.57878    0.63504   0.911  0.36208   
## seasonNH.16.17            -15.34240 1913.69624  -0.008  0.99360   
## seasonSH.10                 2.25150    1.21115   1.859  0.06303 . 
## seasonSH.11               -14.77144 1564.96538  -0.009  0.99247   
## seasonSH.12               -16.21786 2757.06089  -0.006  0.99531   
## seasonSH.13                -0.59348    1.45255  -0.409  0.68285   
## seasonSH.14                 2.18591    1.34308   1.628  0.10362   
## seasonSH.15               -17.17799 2260.24118  -0.008  0.99394   
## seasonSH.16                -0.73902    1.16850  -0.632  0.52709   
## seasonSH.17               -16.09051 2957.38259  -0.005  0.99566   
## anyvac1                    -1.47049    0.59485  -2.472  0.01343 * 
## H3.valid1:seasonNH.10.11  -29.25120 3792.76812  -0.008  0.99385   
## H3.valid1:seasonNH.11.12  -31.42871 3267.08693  -0.010  0.99232   
## H3.valid1:seasonNH.12.13    2.00817 3668.72694   0.001  0.99956   
## H3.valid1:seasonNH.13.14  -12.02767 2957.38320  -0.004  0.99676   
## H3.valid1:seasonNH.14.15    1.33833 3366.90928   0.000  0.99968   
## H3.valid1:seasonNH.15.16  -13.36489 2957.38307  -0.005  0.99639   
## H3.valid1:seasonNH.16.17    1.77251 3522.54827   0.001  0.99960   
## H3.valid1:seasonSH.10     -31.57164 5214.78149  -0.006  0.99517   
## H3.valid1:seasonSH.11       0.66341 3345.92737   0.000  0.99984   
## H3.valid1:seasonSH.12       3.04166 4043.20396   0.001  0.99940   
## H3.valid1:seasonSH.13     -12.87437 2957.38324  -0.004  0.99653   
## H3.valid1:seasonSH.14     -16.06667 2957.38302  -0.005  0.99567   
## H3.valid1:seasonSH.15       2.66549 3722.20402   0.001  0.99943   
## H3.valid1:seasonSH.16     -12.50870 2957.38326  -0.004  0.99663   
## H3.valid1:seasonSH.17            NA         NA      NA       NA   
## H3.valid1:anyvac1           1.30719    0.74334   1.759  0.07865 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 461.40  on 1228  degrees of freedom
## Residual deviance: 372.63  on 1185  degrees of freedom
##   (56 observations deleted due to missingness)
## AIC: 460.63
## 
## Number of Fisher Scoring iterations: 17
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(dth ~ ns(age,knots = 65) + anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit1)
```

```
## 
## Call:
## glm(formula = dth ~ ns(age, knots = 65) + anyav + anydx + country + 
##     H3.valid * season + H3.valid * anyvac, family = binomial, 
##     data = train)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -1.55821  -0.29314  -0.17418  -0.06624   3.06609  
## 
## Coefficients: (1 not defined because of singularities)
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                -6.15511    1.80433  -3.411 0.000647 ***
## ns(age, knots = 65)1        5.47998    1.83195   2.991 0.002778 ** 
## ns(age, knots = 65)2        2.21479    0.77520   2.857 0.004276 ** 
## anyav1                      0.48134    0.45832   1.050 0.293614    
## anydx1                      2.55062    1.03750   2.458 0.013955 *  
## countryAustralia           -2.20051    0.96217  -2.287 0.022194 *  
## countryDenmark             -0.72415    1.16053  -0.624 0.532641    
## countryGermany             -0.17445    1.23660  -0.141 0.887811    
## countryGreece              -1.52670    1.18042  -1.293 0.195891    
## countryOther               -0.60912    1.05111  -0.580 0.562251    
## countrySpain               -2.40187    1.51149  -1.589 0.112045    
## countryThailand            -3.13509    1.15062  -2.725 0.006436 ** 
## countryUK                  -1.25340    1.12524  -1.114 0.265321    
## countryUSA                 -1.90051    1.13099  -1.680 0.092881 .  
## H3.valid1                  13.42262 2919.24624   0.005 0.996331    
## seasonNH.10.11             -0.14424    0.67874  -0.213 0.831705    
## seasonNH.11.12              2.75454    1.63753   1.682 0.092544 .  
## seasonNH.12.13            -15.63302 2139.03829  -0.007 0.994169    
## seasonNH.13.14             -1.20496    1.16995  -1.030 0.303045    
## seasonNH.14.15            -15.51306 1602.73145  -0.010 0.992277    
## seasonNH.15.16              0.06764    0.64919   0.104 0.917015    
## seasonNH.16.17            -15.28922 1905.87232  -0.008 0.993599    
## seasonSH.10                 2.87947    1.29287   2.227 0.025935 *  
## seasonSH.11               -14.44201 1552.47763  -0.009 0.992578    
## seasonSH.12               -14.98445 2827.13201  -0.005 0.995771    
## seasonSH.13                 0.11123    1.55478   0.072 0.942965    
## seasonSH.14                 1.93333    1.39232   1.389 0.164965    
## seasonSH.15               -17.08638 2303.32609  -0.007 0.994081    
## seasonSH.16                -0.86111    1.19815  -0.719 0.472325    
## seasonSH.17               -16.62048 2919.24600  -0.006 0.995457    
## anyvac1                    -1.71022    0.60092  -2.846 0.004427 ** 
## H3.valid1:seasonNH.10.11  -29.64416 3737.61897  -0.008 0.993672    
## H3.valid1:seasonNH.11.12  -32.35849 3206.97191  -0.010 0.991949    
## H3.valid1:seasonNH.12.13    1.52320 3619.04456   0.000 0.999664    
## H3.valid1:seasonNH.13.14  -12.80189 2919.24663  -0.004 0.996501    
## H3.valid1:seasonNH.14.15    0.61185 3330.27728   0.000 0.999853    
## H3.valid1:seasonNH.15.16  -14.14470 2919.24650  -0.005 0.996134    
## H3.valid1:seasonNH.16.17    0.75633 3486.30863   0.000 0.999827    
## H3.valid1:seasonSH.10     -33.66135 5271.88372  -0.006 0.994905    
## H3.valid1:seasonSH.11       0.15525 3306.38559   0.000 0.999963    
## H3.valid1:seasonSH.12       0.84684 4063.82511   0.000 0.999834    
## H3.valid1:seasonSH.13     -14.11499 2919.24672  -0.005 0.996142    
## H3.valid1:seasonSH.14     -16.35868 2919.24646  -0.006 0.995529    
## H3.valid1:seasonSH.15       1.95802 3718.50900   0.001 0.999580    
## H3.valid1:seasonSH.16     -13.23605 2919.24666  -0.005 0.996382    
## H3.valid1:seasonSH.17            NA         NA      NA       NA    
## H3.valid1:anyvac1           1.31580    0.74749   1.760 0.078359 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 461.40  on 1228  degrees of freedom
## Residual deviance: 354.66  on 1183  degrees of freedom
##   (56 observations deleted due to missingness)
## AIC: 446.66
## 
## Number of Fisher Scoring iterations: 17
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(dth ~ ns(age, knots = 65) + p.protected + anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit2)
```

```
## 
## Call:
## glm(formula = dth ~ ns(age, knots = 65) + p.protected + anyav + 
##     anydx + country + H3.valid * season + H3.valid * anyvac, 
##     family = binomial, data = train)
## 
## Deviance Residuals: 
##      Min        1Q    Median        3Q       Max  
## -1.59515  -0.29611  -0.16782  -0.05161   3.09852  
## 
## Coefficients: (1 not defined because of singularities)
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                -6.30028    1.83555  -3.432 0.000598 ***
## ns(age, knots = 65)1        7.45134    2.64708   2.815 0.004879 ** 
## ns(age, knots = 65)2        2.49355    0.85163   2.928 0.003412 ** 
## p.protected                -1.03161    0.76111  -1.355 0.175291    
## anyav1                      0.44852    0.45900   0.977 0.328490    
## anydx1                      2.59133    1.04156   2.488 0.012849 *  
## countryAustralia           -2.22292    0.96880  -2.295 0.021762 *  
## countryDenmark             -0.68923    1.15977  -0.594 0.552326    
## countryGermany             -0.19453    1.23409  -0.158 0.874747    
## countryGreece              -1.42036    1.18020  -1.203 0.228786    
## countryOther               -0.54814    1.05706  -0.519 0.604073    
## countrySpain               -2.40546    1.50915  -1.594 0.110954    
## countryThailand            -3.08445    1.14124  -2.703 0.006878 ** 
## countryUK                  -1.18618    1.12498  -1.054 0.291700    
## countryUSA                 -1.81846    1.12859  -1.611 0.107120    
## H3.valid1                  12.50388 2924.82834   0.004 0.996589    
## seasonNH.10.11             -0.16218    0.67590  -0.240 0.810368    
## seasonNH.11.12              2.81600    1.67011   1.686 0.091773 .  
## seasonNH.12.13            -15.59979 2142.56560  -0.007 0.994191    
## seasonNH.13.14             -1.33697    1.17288  -1.140 0.254327    
## seasonNH.14.15            -15.50638 1611.98539  -0.010 0.992325    
## seasonNH.15.16              0.01075    0.64750   0.017 0.986758    
## seasonNH.16.17            -15.25988 1921.91269  -0.008 0.993665    
## seasonSH.10                 2.88857    1.29016   2.239 0.025160 *  
## seasonSH.11               -14.53157 1549.10059  -0.009 0.992515    
## seasonSH.12               -14.88271 2822.08101  -0.005 0.995792    
## seasonSH.13                 0.01210    1.56633   0.008 0.993837    
## seasonSH.14                 1.90958    1.38939   1.374 0.169316    
## seasonSH.15               -17.10297 2303.14897  -0.007 0.994075    
## seasonSH.16                -0.82105    1.19765  -0.686 0.492998    
## seasonSH.17               -16.59892 2924.82802  -0.006 0.995472    
## anyvac1                    -1.68873    0.59899  -2.819 0.004813 ** 
## H3.valid1:seasonNH.10.11  -29.58298 3696.54186  -0.008 0.993615    
## H3.valid1:seasonNH.11.12  -32.45655 3203.11504  -0.010 0.991915    
## H3.valid1:seasonNH.12.13    1.43759 3625.63206   0.000 0.999684    
## H3.valid1:seasonNH.13.14  -12.68465 2924.82865  -0.004 0.996540    
## H3.valid1:seasonNH.14.15    0.57172 3339.62832   0.000 0.999863    
## H3.valid1:seasonNH.15.16  -14.13859 2924.82852  -0.005 0.996143    
## H3.valid1:seasonNH.16.17    0.70714 3499.76695   0.000 0.999839    
## H3.valid1:seasonSH.10     -33.75359 5258.01595  -0.006 0.994878    
## H3.valid1:seasonSH.11       0.40986 3309.73313   0.000 0.999901    
## H3.valid1:seasonSH.12       0.72095 4064.32800   0.000 0.999858    
## H3.valid1:seasonSH.13     -13.96834 2924.82875  -0.005 0.996189    
## H3.valid1:seasonSH.14     -16.30801 2924.82848  -0.006 0.995551    
## H3.valid1:seasonSH.15       1.97120 3722.78323   0.001 0.999578    
## H3.valid1:seasonSH.16     -13.23113 2924.82867  -0.005 0.996391    
## H3.valid1:seasonSH.17            NA         NA      NA       NA    
## H3.valid1:anyvac1           1.27948    0.74815   1.710 0.087231 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 461.40  on 1228  degrees of freedom
## Residual deviance: 352.62  on 1182  degrees of freedom
##   (56 observations deleted due to missingness)
## AIC: 446.62
## 
## Number of Fisher Scoring iterations: 17
```

```r
par(mfrow = c(1, 2), las = 3, mar = c(6, 6, 6, 3), cex.axis = .7)
library(gam)
plot.gam(fit0, se = T, col = 'darkslateblue', main = 'H1N1 baseline')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-3.png)

```r
plot.gam(fit1, se = T, col = 'darkslateblue', main = 'H1N1 baseline + age')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-4.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-5.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-6.png)

```r
plot.gam(fit2, se = T, col = 'darkslateblue', main = 'H1N1 baseline + age + impriting')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-7.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-8.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-9.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-10.png)

```r
pdf('003Death_constrained.pdf')
plot.gam(fit2, se = T, col = 'darkslateblue', main = 'H1N1 baseline + age + impriting')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

```r
dev.off()
```

```
## RStudioGD 
##         2
```

```r
## Test performance on test data
# Predict probs for each patient
probs0 = predict(fit0, newdata = test, type = 'response')
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-11.png)

```r
probs1 = predict(fit1, newdata = test, type = 'response')
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```r
probs2 = predict(fit2, newdata = test, type = 'response')
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```r
# Create a predict object
pr0 = prediction(probs0, test$dth)
pr1 = prediction(probs1, test$dth)
pr2 = prediction(probs2, test$dth)

# Assess tp and fp
prf0 = performance(pr0, measure = 'tpr', x.measure = 'fpr')
prf1 = performance(pr1, measure = 'tpr', x.measure = 'fpr')
prf2 = performance(pr2, measure = 'tpr', x.measure = 'fpr')

# Calculate AUCs
auc = sapply(list(pr0, pr1, pr2), FUN = function(xx){aa = performance(xx, measure = 'auc')
                                                                      round(aa@y.values[[1]], 3) })

# Plot AUC
par(mfrow = c(1,1))
plot(prf0, lwd = 2, col = cols[3])
plot(prf1, col = cols[2], add = TRUE, lwd = 2)
plot(prf2, col = cols[1], add = TRUE, lwd = 2)
abline(a = 0, b = 1, lty = 2)
legend('bottomright', paste(c('baseline                  AUC = ', 'baseline+age          AUC = ', 'baseline+age+imp  AUC = '), auc), col = cols[3:1], lty = 1)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-12.png)

```r
pdf('003Death_AUC_constrained.pdf')
par(lwd = 3.5, cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
# Plot AUC
par(mfrow = c(1,1))
plot(prf0, col = cols[3])
plot(prf1, col = cols[2], add = TRUE)
plot(prf2, col = cols[1], add = TRUE)
abline(a = 0, b = 1)
legend('bottomright', paste(c('baseline                  AUC = ', 'baseline+age          AUC = ', 'baseline+age+imp  AUC = '), auc), col = cols[3:1], lty = 1, bty = 'n', cex = 1.5)
dev.off()
```

```
## RStudioGD 
##         2
```

```r
AIC = c(fit0$aic, fit1$aic, fit2$aic); names(AIC) =  c('baseline', 'baeline+age', 'baseline+age+imp')
del.AIC = sort(AIC - min(AIC))
del.AIC
```

```
## baseline+age+imp      baeline+age         baseline 
##       0.00000000       0.03697211      14.00784714
```

```r
load('master.AIC.table.constrained.RData')
master.AIC.table['003ICU',  c('baseline', 'baseline+age', 'baseline+age+imp')] = AIC - min(AIC)
save(master.AIC.table, file = 'master.AIC.table.constrained.RData')
rm(master.AIC.table)
```


