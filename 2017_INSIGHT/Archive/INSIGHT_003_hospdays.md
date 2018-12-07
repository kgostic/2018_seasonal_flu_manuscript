---
title: "003 hospdays"
author: "Katie Gostic"
date: "11/1/2017"
output: html_document
---

#### Note: hospdays records the number of days hospitalized after study enrollment. In another script, I test tothosp (the total number of hospitalized days). I try to control for hospitalization unrelated to infelunza by excluding patients with underlying conditions from the tothosp analysis


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




#### Visualize the relationship between age and duration of hospitalization
![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

```
## 
## 	Spearman's rank correlation rho
## 
## data:  valid$age and valid$hospdays
## S = 56530000, p-value = 2.296e-10
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##      rho 
## 0.227336
```

There is a positive correlation between age and duration of hospitalization.


#### Visualize the relationship between age and prob of H3N2 infection
![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

```
## 
## 	Spearman's rank correlation rho
## 
## data:  valid$age and valid$hospdays
## S = 59538000, p-value = 1.403e-07
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.1894313
```

Results show a positive correlation between duration of hosp and age.

#### Visualize the relationship between country and symptom duration
H1N1

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

H3N2

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

# No clear effect of country

#### Visualize the relationship between season and symptom duration

H1N1

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

No clear effect of season.



H3N2

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)


No clear effect of season
--> Try models that don't include season and country as blocking vars


### H1N1 10-fold CV
We know the relationship between age and probability of infection should be non-linear. Use cross-validation to verify, and to choose the number of degrees of freedom to include in a spline


Plot the overall test error rate vs. df


No clear preference. Use 4 df?


### H3N2 10-fold CV
We know the relationship between age and probability of infection should be non-linear. Use cross-validation to verify, and to choose the number of degrees of freedom to include in a spline



Plot the overall test error rate vs. df


Conclusion: Clear preference for linear fit for H3N2

## Set up a logistic regression that includes the following predictors:
* age
* imprinting status
* vaccination status
* antiviral use
* season
* country





### Clean data for fitting

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## 
## Argentina Australia   Austria   Belgium     Chile     China   Denmark 
##        50        69         3        11         1        15        67 
##   Germany    Greece    Norway      Peru    Poland Singapore     Spain 
##        41        81         7        15         6         6        40 
##  Thailand        UK       USA 
##        64       131       114
```

```
## [1] 721
```


### First fit a model where the only independent variable is age
Treat vaccination, antivial use, underlying symptoms, country and season as blocking variables with fixed effects

```r
library(splines)
## H1N1 fit to medical history only (no effect of country, season, age or imprinting)
fit00 = glm(hospdays ~ anyvac + anyav + anydx, data = train.H1N1)
summary(fit00)
```

```
## 
## Call:
## glm(formula = hospdays ~ anyvac + anyav + anydx, data = train.H1N1)
## 
## Deviance Residuals: 
##    Min      1Q  Median      3Q     Max  
## -9.113  -6.113  -3.285   0.715  90.122  
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    5.051      1.264   3.995 7.23e-05 ***
## anyvac1       -2.893      1.028  -2.814  0.00504 ** 
## anyav1         1.235      1.184   1.043  0.29740    
## anydx1         2.828      1.032   2.740  0.00631 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 126.5023)
## 
##     Null deviance: 83312  on 648  degrees of freedom
## Residual deviance: 81594  on 645  degrees of freedom
## AIC: 4989.1
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to country, season and medical history only (no effect of age or imprinting)  
fit0 = glm(hospdays ~ anyvac + anyav + anydx + country + season, data = train.H1N1)
summary(fit0)
```

```
## 
## Call:
## glm(formula = hospdays ~ anyvac + anyav + anydx + country + season, 
##     data = train.H1N1)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -17.162   -5.386   -2.596    1.266   91.161  
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        5.1952     3.1359   1.657 0.098086 .  
## anyvac1           -2.4080     1.0689  -2.253 0.024620 *  
## anyav1             1.5469     1.1930   1.297 0.195255    
## anydx1             1.8164     1.0550   1.722 0.085618 .  
## countryAustralia  -1.5888     2.5730  -0.617 0.537149    
## countryBelgium    -6.1290     4.6112  -1.329 0.184290    
## countryDenmark    -4.2891     2.9824  -1.438 0.150905    
## countryGermany     5.7887     3.2194   1.798 0.072654 .  
## countryGreece     -2.5102     2.9260  -0.858 0.391277    
## countryOther       0.9888     2.9829   0.331 0.740399    
## countrySpain      -0.7366     3.2180  -0.229 0.819026    
## countryThailand   -5.8312     2.5775  -2.262 0.024020 *  
## countryUK         -1.3214     2.7698  -0.477 0.633491    
## countryUSA        -2.3068     2.8077  -0.822 0.411636    
## seasonNH.10.11     6.3614     1.6760   3.796 0.000162 ***
## seasonNH.11.12    -2.3167     4.7205  -0.491 0.623761    
## seasonNH.12.13     2.9316     3.2325   0.907 0.364816    
## seasonNH.13.14     2.1346     1.7764   1.202 0.229957    
## seasonNH.14.15    13.0224     3.4878   3.734 0.000206 ***
## seasonNH.15.16     0.7673     1.5281   0.502 0.615750    
## seasonNH.16.17    -3.3473     3.9887  -0.839 0.401688    
## seasonSH.10        6.8966     3.5137   1.963 0.050120 .  
## seasonSH.11        2.2639     3.4961   0.648 0.517516    
## seasonSH.12       -1.8197     5.8075  -0.313 0.754136    
## seasonSH.13        0.3320     3.4332   0.097 0.922994    
## seasonSH.14        5.4597     4.0784   1.339 0.181165    
## seasonSH.15        0.6925     4.4448   0.156 0.876248    
## seasonSH.16        0.8229     2.4167   0.340 0.733600    
## seasonSH.17       -4.6072     6.5399  -0.704 0.481401    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 118.5167)
## 
##     Null deviance: 83312  on 648  degrees of freedom
## Residual deviance: 73480  on 620  degrees of freedom
## AIC: 4971.1
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(hospdays ~ ns(age, knots = 65) + anyvac + anyav + anydx + country + season, data = train.H1N1)
summary(fit1)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + anyvac + anyav + 
##     anydx + country + season, data = train.H1N1)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -18.981   -5.262   -2.697    1.550   90.474  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           0.96371    3.34046   0.288 0.773064    
## ns(age, knots = 65)1 11.34994    3.19497   3.552 0.000411 ***
## ns(age, knots = 65)2  0.43640    2.65990   0.164 0.869732    
## anyvac1              -2.80204    1.08594  -2.580 0.010102 *  
## anyav1                1.55023    1.18300   1.310 0.190541    
## anydx1                0.83728    1.08774   0.770 0.441744    
## countryAustralia     -1.11212    2.57553  -0.432 0.666035    
## countryBelgium       -4.90543    4.60340  -1.066 0.287016    
## countryDenmark       -3.47328    2.98281  -1.164 0.244699    
## countryGermany        7.10683    3.21981   2.207 0.027665 *  
## countryGreece        -2.62952    2.90142  -0.906 0.365138    
## countryOther          1.26712    2.96640   0.427 0.669414    
## countrySpain         -0.57402    3.19949  -0.179 0.857676    
## countryThailand      -5.39797    2.56083  -2.108 0.035442 *  
## countryUK            -0.35799    2.77696  -0.129 0.897468    
## countryUSA           -1.92188    2.79849  -0.687 0.492494    
## seasonNH.10.11        5.79507    1.67067   3.469 0.000559 ***
## seasonNH.11.12       -0.39534    4.71475  -0.084 0.933201    
## seasonNH.12.13        1.73892    3.22383   0.539 0.589807    
## seasonNH.13.14        1.79960    1.76440   1.020 0.308151    
## seasonNH.14.15       12.48132    3.50335   3.563 0.000395 ***
## seasonNH.15.16        0.08157    1.53523   0.053 0.957646    
## seasonNH.16.17       -3.47768    3.95594  -0.879 0.379687    
## seasonSH.10           7.75423    3.49237   2.220 0.026758 *  
## seasonSH.11           2.63480    3.46956   0.759 0.447901    
## seasonSH.12          -0.13154    5.77864  -0.023 0.981847    
## seasonSH.13           1.17584    3.41244   0.345 0.730530    
## seasonSH.14           5.52103    4.05034   1.363 0.173346    
## seasonSH.15          -0.40962    4.41848  -0.093 0.926168    
## seasonSH.16           0.48568    2.40137   0.202 0.839785    
## seasonSH.17          -6.00971    6.49700  -0.925 0.355328    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 116.5131)
## 
##     Null deviance: 83312  on 648  degrees of freedom
## Residual deviance: 72005  on 618  degrees of freedom
## AIC: 4962
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(hospdays ~ ns(age, knots = 65) + p.g1.protection + anyvac + anyav + anydx + country + season, data = train.H1N1)
summary(fit2)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + p.g1.protection + 
##     anyvac + anyav + anydx + country + season, data = train.H1N1)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -18.769   -5.346   -2.647    1.540   89.890  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            0.9197     3.3352   0.276 0.782837    
## ns(age, knots = 65)1   5.7473     4.5488   1.263 0.206898    
## ns(age, knots = 65)2  -2.0765     3.0279  -0.686 0.493095    
## p.g1.protection        3.0402     1.7598   1.728 0.084554 .  
## anyvac1               -2.7918     1.0842  -2.575 0.010257 *  
## anyav1                 1.6714     1.1832   1.413 0.158266    
## anydx1                 0.7456     1.0873   0.686 0.493160    
## countryAustralia      -1.1423     2.5715  -0.444 0.657038    
## countryBelgium        -5.1160     4.5976  -1.113 0.266254    
## countryDenmark        -3.4605     2.9780  -1.162 0.245677    
## countryGermany         7.1517     3.2148   2.225 0.026465 *  
## countryGreece         -2.8169     2.8988  -0.972 0.331564    
## countryOther           1.1483     2.9624   0.388 0.698423    
## countrySpain          -0.4650     3.1950  -0.146 0.884341    
## countryThailand       -5.3057     2.5573  -2.075 0.038424 *  
## countryUK             -0.5569     2.7749  -0.201 0.841013    
## countryUSA            -2.1257     2.7965  -0.760 0.447465    
## seasonNH.10.11         5.8614     1.6684   3.513 0.000475 ***
## seasonNH.11.12        -0.8228     4.7137  -0.175 0.861478    
## seasonNH.12.13         2.0562     3.2239   0.638 0.523841    
## seasonNH.13.14         2.1344     1.7722   1.204 0.228900    
## seasonNH.14.15        12.9754     3.5094   3.697 0.000237 ***
## seasonNH.15.16         0.5348     1.5551   0.344 0.731025    
## seasonNH.16.17        -3.6149     3.9504  -0.915 0.360509    
## seasonSH.10            7.5846     3.4882   2.174 0.030054 *  
## seasonSH.11            2.8153     3.4656   0.812 0.416903    
## seasonSH.12           -0.3855     5.7713  -0.067 0.946765    
## seasonSH.13            1.1511     3.4070   0.338 0.735580    
## seasonSH.14            5.4317     4.0442   1.343 0.179738    
## seasonSH.15           -0.4424     4.4114  -0.100 0.920154    
## seasonSH.16            0.6740     2.4000   0.281 0.778946    
## seasonSH.17           -6.1310     6.4870  -0.945 0.344961    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 116.1401)
## 
##     Null deviance: 83312  on 648  degrees of freedom
## Residual deviance: 71658  on 617  degrees of freedom
## AIC: 4960.8
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions NO COUNTRY AND SEASON
fit3 = glm(hospdays ~ ns(age, knots = 65) + anyvac + anyav + anydx, data = train.H1N1)
summary(fit3)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + anyvac + anyav + 
##     anydx, data = train.H1N1)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -10.140   -5.957   -3.307    0.328   89.376  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)   
## (Intercept)            2.0665     1.7045   1.212  0.22580   
## ns(age, knots = 65)1   8.7340     3.0634   2.851  0.00450 **
## ns(age, knots = 65)2   0.3014     2.6017   0.116  0.90782   
## anyvac1               -3.1277     1.0383  -3.012  0.00269 **
## anyav1                 1.0916     1.1796   0.925  0.35512   
## anydx1                 2.0489     1.0664   1.921  0.05513 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 125.3116)
## 
##     Null deviance: 83312  on 648  degrees of freedom
## Residual deviance: 80575  on 643  degrees of freedom
## AIC: 4984.9
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, NO COUNTRY AND SEASON
fit4 = glm(hospdays ~ ns(age, knots = 65) + p.g1.protection + anyvac + anyav + anydx, data = train.H1N1)
summary(fit4)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + p.g1.protection + 
##     anyvac + anyav + anydx, data = train.H1N1)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -10.596   -5.887   -3.215    0.473   88.813  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)   
## (Intercept)             2.054      1.701   1.207  0.22776   
## ns(age, knots = 65)1    2.998      4.368   0.686  0.49268   
## ns(age, knots = 65)2   -2.239      2.941  -0.761  0.44690   
## p.g1.protection         3.239      1.761   1.839  0.06636 . 
## anyvac1                -3.139      1.036  -3.029  0.00255 **
## anyav1                  1.269      1.181   1.074  0.28310   
## anydx1                  1.911      1.067   1.791  0.07378 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 124.849)
## 
##     Null deviance: 83312  on 648  degrees of freedom
## Residual deviance: 80153  on 642  degrees of freedom
## AIC: 4983.5
## 
## Number of Fisher Scoring iterations: 2
```

```r
## Plot the fits for each model side by side
par(mfrow = c(1, 2), las = 3, mar = c(6, 6, 6, 3), cex.axis = .7)
library(gam)
plot.gam(fit00, se = T, col = 'dodgerblue', main = 'Simplified Baseline')
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

```r
plot.gam(fit0, se = T, col = 'dodgerblue', main = 'Baseline')
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-2.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-3.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-4.png)

```r
plot.gam(fit1, se = T, col = 'dodgerblue', main = 'Baseline + age + imp')
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-5.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-6.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-7.png)

```r
plot.gam(fit2, se = T, col = 'dodgerblue', main = 'Baseline + age + imp')
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-8.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-9.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-10.png)

```r
plot.gam(fit3, se = T, col = 'dodgerblue', main = 'Simple baseline + age + imp')
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-11.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-12.png)

```r
plot.gam(fit4, se = T, col = 'dodgerblue', main = 'Simple baseline + age + imp')
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-13.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-14.png)![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-15.png)

```r
pdf('003Hospdays_H1N1.pdf')
plot.gam(fit4, se = T, col = 'dodgerblue', main = 'Simple baseline + age + imp')
dev.off()
```

```
## RStudioGD 
##         2
```

```r
# Predict probs for each patient
pred00 = predict(fit00, newdata = (test.H1N1))
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-16.png)

```r
pred0 = predict(fit0, newdata = (test.H1N1))
pred1 = predict(fit1, newdata = (test.H1N1))
pred2 = predict(fit2, newdata = (test.H1N1))
pred3 = predict(fit3, newdata = (test.H1N1))
pred4 = predict(fit4, newdata = (test.H1N1))
  
# Estimate the simulated error rate
MSE = numeric(6); names(MSE) = c('simple.baseline', 'baseline', 'baeline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')
MSE[1] = mean((pred00 - test.H1N1$hospdays)^2)
MSE[2] = mean((pred0 - test.H1N1$hospdays)^2)
MSE[3] = mean((pred1 - test.H1N1$hospdays)^2)
MSE[4] = mean((pred2 - test.H1N1$hospdays)^2)
MSE[5] = mean((pred3 - test.H1N1$hospdays)^2)
MSE[6] = mean((pred4 - test.H1N1$hospdays)^2)
sort(MSE)
```

```
##         simple.baseline     simple.baseline+age simple.baseline+age+imp 
##                193.6096                194.4922                195.8550 
##             baeline+age                baseline        baseline+age+imp 
##                212.1769                213.8063                215.4716
```

```r
AIC = c(fit00$aic, fit0$aic, fit1$aic, fit2$aic, fit3$aic, fit4$aic); names(AIC) =  c('simple.baseline', 'baseline', 'baeline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')
del.AIC = sort(AIC - min(AIC))
del.AIC
```

```
##        baseline+age+imp             baeline+age                baseline 
##                0.000000                1.131955               10.294744 
## simple.baseline+age+imp     simple.baseline+age         simple.baseline 
##               22.705885               24.116006               28.269083
```

```r
load('master.AIC.table.H1N1.RData')
master.AIC.table['003Hospdays',  c('simple.baseline', 'baseline', 'baseline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')] = AIC - min(AIC)
save(master.AIC.table, file = 'master.AIC.table.H1N1.RData')
rm(master.AIC.table)
```




### Clean H3N2 data for fitting

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## integer(0)
```

```
## 
## Argentina Australia   Belgium     China   Denmark   Germany    Greece 
##        97       120         3         7        26        19        24 
##      Peru Singapore     Spain  Thailand        UK       USA 
##         9        21        29       120        89       174
```

```
## [1] 738
```



```r
library(splines)
## H3N2 fit to medical history only (no effect of country, season, age or imprinting)
fit00 = glm(hospdays ~ anyvac + anyav + anydx, data = train.H3N2)
summary(fit00)
```

```
## 
## Call:
## glm(formula = hospdays ~ anyvac + anyav + anydx, data = train.H3N2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -8.294   -6.038   -3.491    0.509  139.398  
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)   3.7477     1.4930   2.510  0.01231 * 
## anyvac1       0.4358     0.9457   0.461  0.64506   
## anyav1       -0.2563     1.1713  -0.219  0.82687   
## anydx1        4.1110     1.2858   3.197  0.00145 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 136.4937)
## 
##     Null deviance: 91865  on 664  degrees of freedom
## Residual deviance: 90222  on 661  degrees of freedom
## AIC: 5162.5
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H3N2 fit to country, season and medical history only (no effect of age or imprinting)  
fit0 = glm(hospdays ~ anyvac + anyav + anydx + country + season, data = train.H3N2)
summary(fit0)
```

```
## 
## Call:
## glm(formula = hospdays ~ anyvac + anyav + anydx + country + season, 
##     data = train.H3N2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -19.175   -5.851   -2.982    0.991  135.965  
## 
## Coefficients:
##                  Estimate Std. Error t value Pr(>|t|)   
## (Intercept)        0.2553     4.8238   0.053  0.95781   
## anyvac1            0.4668     1.0937   0.427  0.66966   
## anyav1             0.6037     1.2429   0.486  0.62732   
## anydx1             3.5241     1.3613   2.589  0.00985 **
## countryAustralia  -0.1639     2.1259  -0.077  0.93855   
## countryBelgium    -1.5334     7.2287  -0.212  0.83207   
## countryDenmark     0.5326     3.4346   0.155  0.87681   
## countryGreece      5.9178     3.5288   1.677  0.09404 . 
## countryOther       1.5435     3.1143   0.496  0.62034   
## countrySingapore  -7.0528     3.1406  -2.246  0.02507 * 
## countrySpain       0.6919     3.3750   0.205  0.83764   
## countryThailand   -2.7757     2.3147  -1.199  0.23091   
## countryUK         -0.2934     2.7300  -0.107  0.91444   
## countryUSA        -3.0369     2.6069  -1.165  0.24447   
## seasonNH.11.12     2.0254     4.4970   0.450  0.65258   
## seasonNH.12.13     5.1351     3.9649   1.295  0.19574   
## seasonNH.13.14     3.1866     4.3242   0.737  0.46144   
## seasonNH.14.15     5.3234     3.9181   1.359  0.17473   
## seasonNH.15.16     2.6742     4.7655   0.561  0.57488   
## seasonNH.16.17     2.8982     3.9540   0.733  0.46384   
## seasonSH.10       23.2481     7.8341   2.968  0.00311 **
## seasonSH.11        3.4676     5.6670   0.612  0.54082   
## seasonSH.12        3.0283     4.8021   0.631  0.52852   
## seasonSH.13        3.1680     5.0611   0.626  0.53157   
## seasonSH.14        5.4563     4.5026   1.212  0.22603   
## seasonSH.15        4.1515     4.3902   0.946  0.34469   
## seasonSH.16        4.2717     4.4407   0.962  0.33644   
## seasonSH.17        6.6516     4.5888   1.450  0.14769   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 133.3869)
## 
##     Null deviance: 91865  on 664  degrees of freedom
## Residual deviance: 84967  on 637  degrees of freedom
## AIC: 5170.6
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H3N2 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(hospdays ~ ns(age, knots = 65) + anyvac + anyav + anydx + country + season, data = train.H3N2)
summary(fit1)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + anyvac + anyav + 
##     anydx + country + season, data = train.H3N2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -19.115   -5.627   -2.890    1.067  135.897  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)   
## (Intercept)           -0.8431     4.9131  -0.172  0.86381   
## ns(age, knots = 65)1   3.6094     3.4553   1.045  0.29660   
## ns(age, knots = 65)2   2.7793     2.1074   1.319  0.18770   
## anyvac1                0.1813     1.1071   0.164  0.87000   
## anyav1                 0.6413     1.2439   0.516  0.60633   
## anydx1                 2.7128     1.4807   1.832  0.06741 . 
## countryAustralia       0.1038     2.1375   0.049  0.96130   
## countryBelgium        -1.4438     7.2271  -0.200  0.84172   
## countryDenmark         0.7531     3.4383   0.219  0.82670   
## countryGreece          5.8580     3.5275   1.661  0.09727 . 
## countryOther           1.9283     3.1227   0.618  0.53712   
## countrySingapore      -6.5597     3.1705  -2.069  0.03895 * 
## countrySpain           0.7797     3.3739   0.231  0.81731   
## countryThailand       -2.6808     2.3178  -1.157  0.24786   
## countryUK              0.1981     2.7467   0.072  0.94252   
## countryUSA            -2.4175     2.6327  -0.918  0.35883   
## seasonNH.11.12         2.0298     4.5018   0.451  0.65222   
## seasonNH.12.13         4.9497     3.9731   1.246  0.21329   
## seasonNH.13.14         3.1187     4.3243   0.721  0.47106   
## seasonNH.14.15         5.3411     3.9216   1.362  0.17369   
## seasonNH.15.16         2.4816     4.7640   0.521  0.60262   
## seasonNH.16.17         2.9205     3.9578   0.738  0.46084   
## seasonSH.10           21.7692     7.8818   2.762  0.00591 **
## seasonSH.11            4.0570     5.6746   0.715  0.47491   
## seasonSH.12            2.8711     4.8060   0.597  0.55046   
## seasonSH.13            3.0573     5.0724   0.603  0.54690   
## seasonSH.14            5.5132     4.5015   1.225  0.22113   
## seasonSH.15            4.1593     4.3875   0.948  0.34349   
## seasonSH.16            4.2158     4.4382   0.950  0.34254   
## seasonSH.17            6.6920     4.5903   1.458  0.14538   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 133.2175)
## 
##     Null deviance: 91865  on 664  degrees of freedom
## Residual deviance: 84593  on 635  degrees of freedom
## AIC: 5171.7
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H3N2 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(hospdays ~ ns(age, knots = 65) + p.g2.protection + anyvac + anyav + anydx + country + season, data = train.H3N2)
summary(fit2)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + p.g2.protection + 
##     anyvac + anyav + anydx + country + season, data = train.H3N2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -19.101   -5.623   -2.863    1.212  135.755  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)   
## (Intercept)            1.4704     5.4656   0.269  0.78800   
## ns(age, knots = 65)1  -0.5039     5.4825  -0.092  0.92680   
## ns(age, knots = 65)2   1.4949     2.4915   0.600  0.54872   
## p.g2.protection       -2.4152     2.4993  -0.966  0.33424   
## anyvac1                0.1384     1.1080   0.125  0.90060   
## anyav1                 0.5228     1.2500   0.418  0.67595   
## anydx1                 2.7479     1.4813   1.855  0.06404 . 
## countryAustralia       0.1349     2.1379   0.063  0.94969   
## countryBelgium        -1.0278     7.2403  -0.142  0.88716   
## countryDenmark         0.7547     3.4385   0.219  0.82634   
## countryGreece          5.7889     3.5284   1.641  0.10137   
## countryOther           1.7808     3.1266   0.570  0.56917   
## countrySingapore      -6.7182     3.1749  -2.116  0.03473 * 
## countrySpain           0.7032     3.3750   0.208  0.83501   
## countryThailand       -2.7598     2.3194  -1.190  0.23453   
## countryUK              0.1803     2.7469   0.066  0.94767   
## countryUSA            -2.4039     2.6329  -0.913  0.36158   
## seasonNH.11.12         2.0714     4.5022   0.460  0.64562   
## seasonNH.12.13         5.0692     3.9752   1.275  0.20270   
## seasonNH.13.14         3.2990     4.3286   0.762  0.44626   
## seasonNH.14.15         5.5592     3.9283   1.415  0.15751   
## seasonNH.15.16         2.6531     4.7676   0.556  0.57808   
## seasonNH.16.17         3.1877     3.9677   0.803  0.42203   
## seasonSH.10           22.1228     7.8907   2.804  0.00521 **
## seasonSH.11            4.4213     5.6875   0.777  0.43722   
## seasonSH.12            2.8132     4.8067   0.585  0.55858   
## seasonSH.13            3.2607     5.0770   0.642  0.52094   
## seasonSH.14            5.5030     4.5018   1.222  0.22201   
## seasonSH.15            4.1719     4.3877   0.951  0.34206   
## seasonSH.16            4.3765     4.4416   0.985  0.32483   
## seasonSH.17            6.8764     4.5945   1.497  0.13498   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 133.2314)
## 
##     Null deviance: 91865  on 664  degrees of freedom
## Residual deviance: 84469  on 634  degrees of freedom
## AIC: 5172.7
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H3N2 fit to age only, plus fixed effects of vaccination, av use, underlying conditions NO COUNTRY AND SEASON
fit3 = glm(hospdays ~ ns(age, knots = 65) + anyvac + anyav + anydx, data = train.H3N2)
summary(fit3)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + anyvac + anyav + 
##     anydx, data = train.H3N2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -10.213   -5.823   -3.025    0.539  139.124  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)  
## (Intercept)            2.7231     1.8396   1.480   0.1393  
## ns(age, knots = 65)1   4.3949     3.4057   1.290   0.1973  
## ns(age, knots = 65)2   4.4872     2.0101   2.232   0.0259 *
## anyvac1                0.2186     0.9461   0.231   0.8173  
## anyav1                -0.2144     1.1685  -0.184   0.8545  
## anydx1                 3.1169     1.3953   2.234   0.0258 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 135.534)
## 
##     Null deviance: 91865  on 664  degrees of freedom
## Residual deviance: 89317  on 659  degrees of freedom
## AIC: 5159.8
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H3N2 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, NO COUNTRY AND SEASON
fit4 = glm(hospdays ~ ns(age, knots = 65) + p.g2.protection + anyvac + anyav + anydx, data = train.H3N2)
summary(fit4)
```

```
## 
## Call:
## glm(formula = hospdays ~ ns(age, knots = 65) + p.g2.protection + 
##     anyvac + anyav + anydx, data = train.H3N2)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -10.182   -5.794   -3.123    0.718  139.059  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)  
## (Intercept)            4.5570     3.0504   1.494   0.1357  
## ns(age, knots = 65)1   1.2382     5.3983   0.229   0.8187  
## ns(age, knots = 65)2   3.4722     2.4199   1.435   0.1518  
## p.g2.protection       -1.8486     2.4523  -0.754   0.4512  
## anyvac1                0.2074     0.9465   0.219   0.8266  
## anyav1                -0.2817     1.1723  -0.240   0.8102  
## anydx1                 3.1573     1.3968   2.260   0.0241 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 135.6228)
## 
##     Null deviance: 91865  on 664  degrees of freedom
## Residual deviance: 89240  on 658  degrees of freedom
## AIC: 5161.2
## 
## Number of Fisher Scoring iterations: 2
```

```r
## Plot the fits for each model side by side
par(mfrow = c(1, 2), las = 3, mar = c(6, 6, 6, 3), cex.axis = .7)
library(gam)
plot.gam(fit00, se = T, col = 'firebrick1', main = 'Simplified Baseline')
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

```r
plot.gam(fit0, se = T, col = 'firebrick1', main = 'Baseline')
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-2.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-3.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-4.png)

```r
plot.gam(fit1, se = T, col = 'firebrick1', main = 'Baseline + age + imp')
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-5.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-6.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-7.png)

```r
plot.gam(fit2, se = T, col = 'firebrick1', main = 'Baseline + age + imp')
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-8.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-9.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-10.png)

```r
plot.gam(fit3, se = T, col = 'firebrick1', main = 'Simple baseline + age + imp')
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-11.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-12.png)

```r
plot.gam(fit4, se = T, col = 'firebrick1', main = 'Simple baseline + age + imp')
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-13.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-14.png)![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-15.png)

```r
pdf('003Hospdays_H3N2.pdf')
plot.gam(fit4, se = T, col = 'firebrick1', main = 'Simple baseline + age + imp')
dev.off()
```

```
## RStudioGD 
##         2
```

```r
# Predict probs for each patient
pred00 = predict(fit00, newdata = test.H3N2)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-16.png)

```r
pred0 = predict(fit0, newdata = test.H3N2)
pred1 = predict(fit1, newdata = test.H3N2)
pred2 = predict(fit2, newdata = test.H3N2)
pred3 = predict(fit3, newdata = test.H3N2)
pred4 = predict(fit4, newdata = test.H3N2)
  
# Estimate the simulated error rate
MSE = numeric(6); names(MSE) = c('simple.baseline', 'baseline', 'baeline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')
MSE[1] = mean((pred00 - test.H3N2$hospdays)^2)
MSE[2] = mean((pred0 - test.H3N2$hospdays)^2)
MSE[3] = mean((pred1 - test.H3N2$hospdays)^2)
MSE[4] = mean((pred2 - test.H3N2$hospdays)^2)
MSE[5] = mean((pred3 - test.H3N2$hospdays)^2)
MSE[6] = mean((pred4 - test.H3N2$hospdays)^2)
sort(MSE)
```

```
## simple.baseline+age+imp     simple.baseline+age         simple.baseline 
##                85.99072                86.24900                87.38044 
##        baseline+age+imp             baeline+age                baseline 
##                94.15226                94.29379                96.36646
```

```r
AIC = c(fit00$aic, fit0$aic, fit1$aic, fit2$aic, fit3$aic, fit4$aic); names(AIC) =  c('simple.baseline', 'baseline', 'baseline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')
del.AIC = sort(AIC - min(AIC))
del.AIC
```

```
##     simple.baseline+age simple.baseline+age+imp         simple.baseline 
##                0.000000                1.425935                2.707376 
##                baseline            baseline+age        baseline+age+imp 
##               10.801896               11.865403               12.886637
```

```r
load('master.AIC.table.H3N2.RData')
master.AIC.table['003Hospdays',  c('simple.baseline', 'baseline', 'baseline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')] = AIC - min(AIC)
save(master.AIC.table, file = 'master.AIC.table.H3N2.RData')
rm(master.AIC.table)
```

