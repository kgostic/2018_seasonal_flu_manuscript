---
title: "003 ICU"
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



## Set up a logistic regression that includes the following predictors:
* age
* imprinting status
* vaccination status
* antiviral use
* season
* country





### Clean data for fitting

```
## 
## Argentina Australia   Austria   Belgium     Chile     China   Denmark 
##        51        70         8        11         1        16        71 
##   Germany    Greece    Norway      Peru    Poland Singapore     Spain 
##        44        90         7        17         6         7        43 
##  Thailand        UK       USA 
##        64       145       120
```

```
## [1] 771
```

```
## 
## Argentina Australia   Belgium     China   Denmark   Germany    Greece 
##       117       123         5         7        32        21        24 
##      Peru Singapore     Spain  Thailand        UK       USA 
##         9        24        31       125        93       180
```

```
## [1] 791
```


### First fit a model where the only independent variable is age
Treat vaccination, antivial use, underlying symptoms, country and season as blocking variables with fixed effects

```r
library(splines)
## H1N1 fit to country, season and medical history only (no effect of age or imprinting)
fit0 = glm(anyicu ~ anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit0)
```

```
## 
## Call:
## glm(formula = anyicu ~ anyav + anydx + country + H3.valid * season + 
##     H3.valid * anyvac, family = binomial, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.4503  -0.5751  -0.4391  -0.2577   2.7712  
## 
## Coefficients: (1 not defined because of singularities)
##                           Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                -0.3828     0.6543  -0.585 0.558499    
## anyav1                      0.3367     0.2313   1.455 0.145568    
## anydx1                      0.5924     0.2465   2.404 0.016237 *  
## countryAustralia           -1.8526     0.4509  -4.109 3.98e-05 ***
## countryAustria             16.7893  1692.1913   0.010 0.992084    
## countryChile               16.8015  2399.5449   0.007 0.994413    
## countryChina               -0.3087     0.7983  -0.387 0.698952    
## countryDenmark             -1.9778     0.6069  -3.259 0.001119 ** 
## countryGermany             -0.6756     0.6451  -1.047 0.294980    
## countryGreece              -2.4849     0.6781  -3.665 0.000248 ***
## countryNorway              -1.9840     1.2503  -1.587 0.112534    
## countryOther               -1.4215     0.5572  -2.551 0.010737 *  
## countryPeru                 1.2260     0.8727   1.405 0.160050    
## countryPoland             -16.7259   970.4980  -0.017 0.986250    
## countrySingapore          -14.3930   768.4467  -0.019 0.985057    
## countrySpain               -1.9348     0.6622  -2.922 0.003478 ** 
## countryThailand            -3.0543     0.5793  -5.273 1.34e-07 ***
## countryUK                  -1.7518     0.5437  -3.222 0.001273 ** 
## countryUSA                 -1.4764     0.5330  -2.770 0.005601 ** 
## H3.valid1                  10.7450   963.0670   0.011 0.991098    
## seasonNH.10.11              0.3843     0.4088   0.940 0.347253    
## seasonNH.11.12             -0.6896     1.2186  -0.566 0.571481    
## seasonNH.12.13              0.1564     0.7227   0.216 0.828674    
## seasonNH.13.14             -0.2923     0.4463  -0.655 0.512409    
## seasonNH.14.15              1.3760     0.7484   1.839 0.065969 .  
## seasonNH.15.16             -0.4942     0.4176  -1.183 0.236655    
## seasonNH.16.17            -13.5416   755.1175  -0.018 0.985692    
## seasonSH.10                 0.2662     0.9391   0.283 0.776832    
## seasonSH.11                -0.1923     0.8575  -0.224 0.822598    
## seasonSH.12                -2.5909     1.4599  -1.775 0.075944 .  
## seasonSH.13                -1.5673     1.0307  -1.521 0.128361    
## seasonSH.14                 0.1057     1.0246   0.103 0.917869    
## seasonSH.15                -0.4659     1.0902  -0.427 0.669115    
## seasonSH.16                -0.4420     0.6062  -0.729 0.465872    
## seasonSH.17               -12.6724   963.0669  -0.013 0.989501    
## anyvac1                    -0.5894     0.2695  -2.187 0.028763 *  
## H3.valid1:seasonNH.10.11  -11.3513   963.0674  -0.012 0.990596    
## H3.valid1:seasonNH.11.12  -10.9476   963.0680  -0.011 0.990930    
## H3.valid1:seasonNH.12.13  -10.5284   963.0673  -0.011 0.991278    
## H3.valid1:seasonNH.13.14  -11.4417   963.0673  -0.012 0.990521    
## H3.valid1:seasonNH.14.15  -12.9616   963.0672  -0.013 0.989262    
## H3.valid1:seasonNH.15.16  -10.7366   963.0677  -0.011 0.991105    
## H3.valid1:seasonNH.16.17    1.9827  1223.8056   0.002 0.998707    
## H3.valid1:seasonSH.10     -26.3554  1524.9606  -0.017 0.986211    
## H3.valid1:seasonSH.11     -11.3239   963.0675  -0.012 0.990619    
## H3.valid1:seasonSH.12      -8.5357   963.0682  -0.009 0.992928    
## H3.valid1:seasonSH.13     -10.6215   963.0681  -0.011 0.991200    
## H3.valid1:seasonSH.14     -12.1919   963.0674  -0.013 0.989900    
## H3.valid1:seasonSH.15     -10.9966   963.0675  -0.011 0.990890    
## H3.valid1:seasonSH.16     -11.1265   963.0673  -0.012 0.990782    
## H3.valid1:seasonSH.17           NA         NA      NA       NA    
## H3.valid1:anyvac1           0.3451     0.3718   0.928 0.353331    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1180.5  on 1349  degrees of freedom
## Residual deviance: 1008.1  on 1299  degrees of freedom
##   (56 observations deleted due to missingness)
## AIC: 1110.1
## 
## Number of Fisher Scoring iterations: 15
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(anyicu ~ ns(age,knots = 65) + anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit1)
```

```
## 
## Call:
## glm(formula = anyicu ~ ns(age, knots = 65) + anyav + anydx + 
##     country + H3.valid * season + H3.valid * anyvac, family = binomial, 
##     data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.4913  -0.5743  -0.4263  -0.2521   2.9142  
## 
## Coefficients: (1 not defined because of singularities)
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                -0.80590    0.69246  -1.164 0.244498    
## ns(age, knots = 65)1        1.07715    0.62688   1.718 0.085744 .  
## ns(age, knots = 65)2       -0.60636    0.48579  -1.248 0.211966    
## anyav1                      0.33478    0.23217   1.442 0.149311    
## anydx1                      0.52026    0.25228   2.062 0.039188 *  
## countryAustralia           -1.89788    0.45568  -4.165 3.11e-05 ***
## countryAustria             16.80599 1669.46079   0.010 0.991968    
## countryChile               16.77562 2399.54485   0.007 0.994422    
## countryChina               -0.40235    0.80117  -0.502 0.615529    
## countryDenmark             -1.99624    0.60929  -3.276 0.001052 ** 
## countryGermany             -0.63846    0.64766  -0.986 0.324236    
## countryGreece              -2.59694    0.68045  -3.817 0.000135 ***
## countryNorway              -2.13041    1.25379  -1.699 0.089285 .  
## countryOther               -1.48639    0.55824  -2.663 0.007754 ** 
## countryPeru                 1.17812    0.88058   1.338 0.180929    
## countryPoland             -16.67466  969.19659  -0.017 0.986273    
## countrySingapore          -14.55409  768.98707  -0.019 0.984900    
## countrySpain               -2.01424    0.66340  -3.036 0.002396 ** 
## countryThailand            -3.08314    0.58125  -5.304 1.13e-07 ***
## countryUK                  -1.79211    0.54974  -3.260 0.001114 ** 
## countryUSA                 -1.53791    0.53784  -2.859 0.004244 ** 
## H3.valid1                  10.88320  962.18006   0.011 0.990975    
## seasonNH.10.11              0.33048    0.41034   0.805 0.420597    
## seasonNH.11.12             -0.50102    1.20830  -0.415 0.678399    
## seasonNH.12.13              0.06796    0.72617   0.094 0.925434    
## seasonNH.13.14             -0.32642    0.44922  -0.727 0.467443    
## seasonNH.14.15              1.39090    0.75507   1.842 0.065462 .  
## seasonNH.15.16             -0.55169    0.42019  -1.313 0.189200    
## seasonNH.16.17            -13.59822  751.88289  -0.018 0.985571    
## seasonSH.10                 0.30653    0.94217   0.325 0.744920    
## seasonSH.11                -0.17976    0.85842  -0.209 0.834132    
## seasonSH.12                -2.40984    1.45913  -1.652 0.098624 .  
## seasonSH.13                -1.51662    1.02841  -1.475 0.140287    
## seasonSH.14                 0.04530    1.02625   0.044 0.964794    
## seasonSH.15                -0.64582    1.09263  -0.591 0.554472    
## seasonSH.16                -0.51816    0.60860  -0.851 0.394546    
## seasonSH.17               -12.81861  962.18003  -0.013 0.989371    
## anyvac1                    -0.61530    0.27323  -2.252 0.024328 *  
## H3.valid1:seasonNH.10.11  -11.32614  962.18046  -0.012 0.990608    
## H3.valid1:seasonNH.11.12  -11.28953  962.18108  -0.012 0.990638    
## H3.valid1:seasonNH.12.13  -10.61555  962.18037  -0.011 0.991197    
## H3.valid1:seasonNH.13.14  -11.52837  962.18043  -0.012 0.990440    
## H3.valid1:seasonNH.14.15  -13.14903  962.18034  -0.014 0.989097    
## H3.valid1:seasonNH.15.16  -10.75722  962.18082  -0.011 0.991080    
## H3.valid1:seasonNH.16.17    1.87282 1221.11349   0.002 0.998776    
## H3.valid1:seasonSH.10     -26.33900 1524.62533  -0.017 0.986217    
## H3.valid1:seasonSH.11     -11.49938  962.18057  -0.012 0.990464    
## H3.valid1:seasonSH.12      -8.87443  962.18134  -0.009 0.992641    
## H3.valid1:seasonSH.13     -10.91505  962.18119  -0.011 0.990949    
## H3.valid1:seasonSH.14     -12.26798  962.18055  -0.013 0.989827    
## H3.valid1:seasonSH.15     -10.94208  962.18061  -0.011 0.990927    
## H3.valid1:seasonSH.16     -11.20010  962.18036  -0.012 0.990713    
## H3.valid1:seasonSH.17            NA         NA      NA       NA    
## H3.valid1:anyvac1           0.40776    0.37393   1.090 0.275502    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1180.5  on 1349  degrees of freedom
## Residual deviance: 1003.2  on 1297  degrees of freedom
##   (56 observations deleted due to missingness)
## AIC: 1109.2
## 
## Number of Fisher Scoring iterations: 15
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(anyicu ~ ns(age, knots = 65) + p.protected + anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit2)
```

```
## 
## Call:
## glm(formula = anyicu ~ ns(age, knots = 65) + p.protected + anyav + 
##     anydx + country + H3.valid * season + H3.valid * anyvac, 
##     family = binomial, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.4964  -0.5765  -0.4255  -0.2519   2.8876  
## 
## Coefficients: (1 not defined because of singularities)
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                -0.86111    0.70163  -1.227 0.219712    
## ns(age, knots = 65)1        1.03555    0.62833   1.648 0.099332 .  
## ns(age, knots = 65)2       -0.57526    0.48920  -1.176 0.239624    
## p.protected                 0.12537    0.24317   0.516 0.606152    
## anyav1                      0.34057    0.23258   1.464 0.143108    
## anydx1                      0.51563    0.25238   2.043 0.041046 *  
## countryAustralia           -1.89301    0.45623  -4.149 3.34e-05 ***
## countryAustria             16.81000 1663.73139   0.010 0.991938    
## countryChile               16.86432 2399.54486   0.007 0.994392    
## countryChina               -0.42357    0.80238  -0.528 0.597576    
## countryDenmark             -1.98141    0.60992  -3.249 0.001159 ** 
## countryGermany             -0.62667    0.64811  -0.967 0.333591    
## countryGreece              -2.61814    0.68200  -3.839 0.000124 ***
## countryNorway              -2.13786    1.25504  -1.703 0.088488 .  
## countryOther               -1.47329    0.55861  -2.637 0.008353 ** 
## countryPeru                 1.19368    0.88270   1.352 0.176280    
## countryPoland             -16.65222  967.34838  -0.017 0.986266    
## countrySingapore          -14.57696  768.76442  -0.019 0.984872    
## countrySpain               -2.00916    0.66372  -3.027 0.002469 ** 
## countryThailand            -3.08301    0.58191  -5.298 1.17e-07 ***
## countryUK                  -1.79055    0.55015  -3.255 0.001135 ** 
## countryUSA                 -1.54475    0.53838  -2.869 0.004114 ** 
## H3.valid1                  10.95429  961.05837   0.011 0.990906    
## seasonNH.10.11              0.32187    0.41075   0.784 0.433258    
## seasonNH.11.12             -0.49823    1.20768  -0.413 0.679938    
## seasonNH.12.13              0.05515    0.72840   0.076 0.939643    
## seasonNH.13.14             -0.31906    0.44962  -0.710 0.477938    
## seasonNH.14.15              1.38040    0.75551   1.827 0.067681 .  
## seasonNH.15.16             -0.55140    0.42022  -1.312 0.189459    
## seasonNH.16.17            -13.61825  751.22815  -0.018 0.985537    
## seasonSH.10                 0.30206    0.94344   0.320 0.748838    
## seasonSH.11                -0.17299    0.85891  -0.201 0.840377    
## seasonSH.12                -2.41869    1.45927  -1.657 0.097425 .  
## seasonSH.13                -1.52372    1.02940  -1.480 0.138820    
## seasonSH.14                 0.02945    1.02349   0.029 0.977042    
## seasonSH.15                -0.65841    1.09060  -0.604 0.546034    
## seasonSH.16                -0.53615    0.61031  -0.878 0.379679    
## seasonSH.17               -12.83704  961.05832  -0.013 0.989343    
## anyvac1                    -0.63553    0.27625  -2.301 0.021418 *  
## H3.valid1:seasonNH.10.11  -11.34827  961.05876  -0.012 0.990579    
## H3.valid1:seasonNH.11.12  -11.31614  961.05937  -0.012 0.990605    
## H3.valid1:seasonNH.12.13  -10.62119  961.05866  -0.011 0.991182    
## H3.valid1:seasonNH.13.14  -11.55811  961.05872  -0.012 0.990405    
## H3.valid1:seasonNH.14.15  -13.16510  961.05863  -0.014 0.989070    
## H3.valid1:seasonNH.15.16  -10.78332  961.05911  -0.011 0.991048    
## H3.valid1:seasonNH.16.17    1.85886 1219.82650   0.002 0.998784    
## H3.valid1:seasonSH.10     -26.36770 1525.44499  -0.017 0.986209    
## H3.valid1:seasonSH.11     -11.56489  961.05887  -0.012 0.990399    
## H3.valid1:seasonSH.12      -8.87195  961.05963  -0.009 0.992634    
## H3.valid1:seasonSH.13     -10.90917  961.05948  -0.011 0.990943    
## H3.valid1:seasonSH.14     -12.26710  961.05884  -0.013 0.989816    
## H3.valid1:seasonSH.15     -10.94467  961.05890  -0.011 0.990914    
## H3.valid1:seasonSH.16     -11.21060  961.05865  -0.012 0.990693    
## H3.valid1:seasonSH.17            NA         NA      NA       NA    
## H3.valid1:anyvac1           0.43818    0.37848   1.158 0.246971    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 1180.5  on 1349  degrees of freedom
## Residual deviance: 1002.9  on 1296  degrees of freedom
##   (56 observations deleted due to missingness)
## AIC: 1110.9
## 
## Number of Fisher Scoring iterations: 15
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
pdf('003ICU_constrained.pdf')
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
pr0 = prediction(probs0, test$anyicu)
pr1 = prediction(probs1, test$anyicu)
pr2 = prediction(probs2, test$anyicu)

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
pdf('003ICU_AUC_constrained.pdf')
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
##      baeline+age         baseline baseline+age+imp 
##        0.0000000        0.8619876        1.7346746
```

```r
load('master.AIC.table.constrained.RData')
master.AIC.table['003ICU',  c('baseline', 'baseline+age', 'baseline+age+imp')] = AIC - min(AIC)
save(master.AIC.table, file = 'master.AIC.table.constrained.RData')
rm(master.AIC.table)
```

