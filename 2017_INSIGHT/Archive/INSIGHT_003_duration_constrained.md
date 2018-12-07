---
title: "003_symdur"
author: "Katie Gostic"
date: "11/1/2017"
output: html_document
---

#### Note: symdur only contains info on the duration of symptoms prior to admission, so all these analyses may be useless. Seems like much of the total duration of symptoms data is missing.


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


Notes on duration of symptoms:
Currently using symdur here, because I want to do this analysis on the continuous duration of symptoms.
But I'm waiting for an email from Deb clarifying whether or not this variable represents total duration of symptoms or just duration prior to enrollment.



#### Visualize the relationship between age and duration of symptoms
![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

```
## 
## 	Spearman's rank correlation rho
## 
## data:  valid$age and valid$symdur
## S = 66961000, p-value = 0.05384
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## 0.07015608
```

There is a trend toward a very weak positive correlation here, but looking at the data shows that it's not very strong.


#### Visualize the relationship between age and prob of H3N2 infection
![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

```
## 
## 	Spearman's rank correlation rho
## 
## data:  valid$age and valid$symdur
## S = 63158000, p-value = 1.29e-05
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.1568748
```

Results suggest a weak, positive correlation, although this may be influenced by outliers.

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
## H1N1 fit to medical history only (no effect of country, season, age or imprinting)
fit00 = glm(symdur ~ anyav + anydx  + H3.valid*anyvac, data = train)
summary(fit00)
```

```
## 
## Call:
## glm(formula = symdur ~ anyav + anydx + H3.valid * anyvac, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -7.577   -3.063   -1.255    0.839  143.647  
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         7.1610     0.5883  12.173  < 2e-16 ***
## anyav1             -1.5140     0.4947  -3.060  0.00226 ** 
## anydx1              0.4163     0.4626   0.900  0.36832    
## H3.valid1          -0.8081     0.4619  -1.749  0.08046 .  
## anyvac1             0.3356     0.5998   0.560  0.57586    
## H3.valid1:anyvac1   0.6472     0.7834   0.826  0.40888    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 44.25767)
## 
##     Null deviance: 58890  on 1320  degrees of freedom
## Residual deviance: 58199  on 1315  degrees of freedom
##   (85 observations deleted due to missingness)
## AIC: 8763.4
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to country, season and medical history only (no effect of age or imprinting)  
fit0 = glm(symdur ~ anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train)
summary(fit0)
```

```
## 
## Call:
## glm(formula = symdur ~ anyav + anydx + country + H3.valid * season + 
##     H3.valid * anyvac, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -8.530   -2.942   -1.183    1.130  139.763  
## 
## Coefficients: (1 not defined because of singularities)
##                          Estimate Std. Error t value Pr(>|t|)    
## (Intercept)               4.91286    1.35989   3.613 0.000315 ***
## anyav1                   -1.00594    0.51824  -1.941 0.052473 .  
## anydx1                    0.16581    0.48272   0.343 0.731295    
## countryAustralia          0.85222    0.95004   0.897 0.369871    
## countryBelgium            2.70497    2.63013   1.028 0.303931    
## countryDenmark            3.47063    1.27615   2.720 0.006624 ** 
## countryGermany            2.10743    1.41954   1.485 0.137899    
## countryGreece             0.32173    1.24783   0.258 0.796577    
## countryOther             -0.25165    1.30681  -0.193 0.847329    
## countrySpain              1.51848    1.33407   1.138 0.255235    
## countryThailand          -1.76516    0.97379  -1.813 0.070117 .  
## countryUK                 1.66027    1.12775   1.472 0.141214    
## countryUSA                0.21228    1.11773   0.190 0.849402    
## countrySingapore         -0.96518    1.67268  -0.577 0.564025    
## H3.valid1                -1.06552    3.20459  -0.332 0.739566    
## seasonNH.10.11            1.34408    0.99993   1.344 0.179130    
## seasonNH.11.12            0.09015    2.82800   0.032 0.974576    
## seasonNH.12.13            0.40384    2.15956   0.187 0.851689    
## seasonNH.13.14            1.59195    0.98443   1.617 0.106097    
## seasonNH.14.15            1.80856    2.05318   0.881 0.378559    
## seasonNH.15.16            1.07776    0.87894   1.226 0.220348    
## seasonNH.16.17            1.04766    2.47814   0.423 0.672541    
## seasonSH.10               1.56159    1.91555   0.815 0.415100    
## seasonSH.11               1.97540    1.89386   1.043 0.297121    
## seasonSH.12               1.71241    3.11734   0.549 0.582883    
## seasonSH.13               3.22883    1.83528   1.759 0.078765 .  
## seasonSH.14               5.32203    2.35334   2.261 0.023898 *  
## seasonSH.15               0.83039    2.73038   0.304 0.761079    
## seasonSH.16               1.70888    1.22275   1.398 0.162485    
## seasonSH.17               2.41749    3.09386   0.781 0.434722    
## anyvac1                   0.13171    0.63001   0.209 0.834436    
## H3.valid1:seasonNH.10.11 -0.60816    4.08288  -0.149 0.881613    
## H3.valid1:seasonNH.11.12  3.04085    4.48222   0.678 0.497626    
## H3.valid1:seasonNH.12.13  3.25205    3.89550   0.835 0.403975    
## H3.valid1:seasonNH.13.14  0.27736    3.48000   0.080 0.936488    
## H3.valid1:seasonNH.14.15  1.11051    3.78280   0.294 0.769136    
## H3.valid1:seasonNH.15.16  1.36841    3.56393   0.384 0.701072    
## H3.valid1:seasonNH.16.17  0.93201    4.01873   0.232 0.816640    
## H3.valid1:seasonSH.10    -1.01020    5.05185  -0.200 0.841538    
## H3.valid1:seasonSH.11    -1.21035    4.06798  -0.298 0.766109    
## H3.valid1:seasonSH.12     0.93840    4.55941   0.206 0.836967    
## H3.valid1:seasonSH.13    -2.72823    4.10624  -0.664 0.506548    
## H3.valid1:seasonSH.14    -2.62183    3.96045  -0.662 0.508090    
## H3.valid1:seasonSH.15     1.43524    4.19191   0.342 0.732118    
## H3.valid1:seasonSH.16     0.53177    3.45927   0.154 0.877852    
## H3.valid1:seasonSH.17          NA         NA      NA       NA    
## H3.valid1:anyvac1        -0.21749    0.82925  -0.262 0.793156    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 43.82223)
## 
##     Null deviance: 58890  on 1320  degrees of freedom
## Residual deviance: 55873  on 1275  degrees of freedom
##   (85 observations deleted due to missingness)
## AIC: 8789.6
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(symdur ~ ns(age, knots = 65) + anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train)
summary(fit1)
```

```
## 
## Call:
## glm(formula = symdur ~ ns(age, knots = 65) + anyav + anydx + 
##     country + H3.valid * season + H3.valid * anyvac, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -8.851   -2.888   -1.138    1.161  140.236  
## 
## Coefficients: (1 not defined because of singularities)
##                          Estimate Std. Error t value Pr(>|t|)   
## (Intercept)               4.27185    1.43402   2.979  0.00295 **
## ns(age, knots = 65)1      1.88012    1.33824   1.405  0.16029   
## ns(age, knots = 65)2     -0.18533    1.03234  -0.180  0.85755   
## anyav1                   -1.02046    0.51849  -1.968  0.04927 * 
## anydx1                   -0.03047    0.50538  -0.060  0.95193   
## countryAustralia          0.86428    0.95326   0.907  0.36476   
## countryBelgium            2.86291    2.63548   1.086  0.27755   
## countryDenmark            3.53969    1.28079   2.764  0.00580 **
## countryGermany            2.24812    1.42451   1.578  0.11478   
## countryGreece             0.25382    1.24871   0.203  0.83896   
## countryOther             -0.22126    1.30723  -0.169  0.86562   
## countrySpain              1.50628    1.33463   1.129  0.25927   
## countryThailand          -1.74931    0.97451  -1.795  0.07288 . 
## countryUK                 1.74768    1.13439   1.541  0.12366   
## countryUSA                0.27961    1.12659   0.248  0.80403   
## countrySingapore         -1.07535    1.68567  -0.638  0.52363   
## H3.valid1                -0.88736    3.21261  -0.276  0.78243   
## seasonNH.10.11            1.26096    1.00180   1.259  0.20837   
## seasonNH.11.12            0.23923    2.83063   0.085  0.93266   
## seasonNH.12.13            0.18945    2.16518   0.087  0.93029   
## seasonNH.13.14            1.51546    0.98700   1.535  0.12493   
## seasonNH.14.15            1.79118    2.05771   0.870  0.38421   
## seasonNH.15.16            0.95122    0.88586   1.074  0.28313   
## seasonNH.16.17            0.95503    2.47920   0.385  0.70014   
## seasonSH.10               1.73497    1.91941   0.904  0.36622   
## seasonSH.11               2.03944    1.89441   1.077  0.28188   
## seasonSH.12               1.98314    3.12311   0.635  0.52555   
## seasonSH.13               3.32916    1.83660   1.813  0.07012 . 
## seasonSH.14               5.31701    2.35494   2.258  0.02413 * 
## seasonSH.15               0.64971    2.73368   0.238  0.81218   
## seasonSH.16               1.62276    1.22503   1.325  0.18552   
## seasonSH.17               2.21666    3.09707   0.716  0.47429   
## anyvac1                   0.06063    0.63529   0.095  0.92399   
## H3.valid1:seasonNH.10.11 -0.72424    4.08364  -0.177  0.85926   
## H3.valid1:seasonNH.11.12  2.64449    4.49214   0.589  0.55617   
## H3.valid1:seasonNH.12.13  3.16143    3.89647   0.811  0.41731   
## H3.valid1:seasonNH.13.14  0.09453    3.48283   0.027  0.97835   
## H3.valid1:seasonNH.14.15  0.86916    3.79178   0.229  0.81873   
## H3.valid1:seasonNH.15.16  1.32496    3.56404   0.372  0.71013   
## H3.valid1:seasonNH.16.17  0.76080    4.02161   0.189  0.84998   
## H3.valid1:seasonSH.10    -1.38773    5.06451  -0.274  0.78412   
## H3.valid1:seasonSH.11    -1.45842    4.07303  -0.358  0.72035   
## H3.valid1:seasonSH.12     0.36711    4.57690   0.080  0.93608   
## H3.valid1:seasonSH.13    -3.22982    4.12299  -0.783  0.43356   
## H3.valid1:seasonSH.14    -2.80036    3.96434  -0.706  0.48008   
## H3.valid1:seasonSH.15     1.44246    4.19216   0.344  0.73084   
## H3.valid1:seasonSH.16     0.39476    3.46079   0.114  0.90920   
## H3.valid1:seasonSH.17          NA         NA      NA       NA   
## H3.valid1:anyvac1        -0.19416    0.82979  -0.234  0.81503   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 43.82091)
## 
##     Null deviance: 58890  on 1320  degrees of freedom
## Residual deviance: 55784  on 1273  degrees of freedom
##   (85 observations deleted due to missingness)
## AIC: 8791.5
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(symdur ~ ns(age, knots = 65) +p.protected + anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train)
summary(fit2)
```

```
## 
## Call:
## glm(formula = symdur ~ ns(age, knots = 65) + p.protected + anyav + 
##     anydx + country + H3.valid * season + H3.valid * anyvac, 
##     data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -8.777   -2.888   -1.133    1.187  140.358  
## 
## Coefficients: (1 not defined because of singularities)
##                          Estimate Std. Error t value Pr(>|t|)   
## (Intercept)               4.37357    1.46072   2.994  0.00281 **
## ns(age, knots = 65)1      1.89054    1.33899   1.412  0.15822   
## ns(age, knots = 65)2     -0.25080    1.04780  -0.239  0.81086   
## p.protected              -0.19497    0.52803  -0.369  0.71201   
## anyav1                   -1.03057    0.51939  -1.984  0.04745 * 
## anydx1                   -0.02501    0.50577  -0.049  0.96058   
## countryAustralia          0.85682    0.95380   0.898  0.36918   
## countryBelgium            2.82224    2.63868   1.070  0.28502   
## countryDenmark            3.52471    1.28187   2.750  0.00605 **
## countryGermany            2.22714    1.42613   1.562  0.11862   
## countryGreece             0.27331    1.25025   0.219  0.82700   
## countryOther             -0.21617    1.30775  -0.165  0.86874   
## countrySpain              1.49393    1.33550   1.119  0.26351   
## countryThailand          -1.75786    0.97512  -1.803  0.07167 . 
## countryUK                 1.74598    1.13478   1.539  0.12415   
## countryUSA                0.29374    1.12763   0.260  0.79453   
## countrySingapore         -1.10724    1.68845  -0.656  0.51209   
## H3.valid1                -0.97313    3.22208  -0.302  0.76269   
## seasonNH.10.11            1.27720    1.00311   1.273  0.20316   
## seasonNH.11.12            0.23374    2.83163   0.083  0.93423   
## seasonNH.12.13            0.22009    2.16750   0.102  0.91914   
## seasonNH.13.14            1.51046    0.98743   1.530  0.12634   
## seasonNH.14.15            1.81417    2.05935   0.881  0.37852   
## seasonNH.15.16            0.95985    0.88647   1.083  0.27911   
## seasonNH.16.17            0.98544    2.48141   0.397  0.69134   
## seasonSH.10               1.72699    1.92019   0.899  0.36862   
## seasonSH.11               2.02484    1.89547   1.068  0.28561   
## seasonSH.12               1.96213    3.12469   0.628  0.53015   
## seasonSH.13               3.31939    1.83742   1.807  0.07107 . 
## seasonSH.14               5.34294    2.35678   2.267  0.02355 * 
## seasonSH.15               0.69104    2.73690   0.252  0.80070   
## seasonSH.16               1.64535    1.22697   1.341  0.18016   
## seasonSH.17               2.24989    3.09943   0.726  0.46803   
## anyvac1                   0.09492    0.64226   0.148  0.88253   
## H3.valid1:seasonNH.10.11 -0.72465    4.08503  -0.177  0.85923   
## H3.valid1:seasonNH.11.12  2.67496    4.49442   0.595  0.55183   
## H3.valid1:seasonNH.12.13  3.15079    3.89790   0.808  0.41905   
## H3.valid1:seasonNH.13.14  0.12688    3.48511   0.036  0.97096   
## H3.valid1:seasonNH.14.15  0.88255    3.79324   0.233  0.81606   
## H3.valid1:seasonNH.15.16  1.35969    3.56649   0.381  0.70309   
## H3.valid1:seasonNH.16.17  0.77555    4.02317   0.193  0.84717   
## H3.valid1:seasonSH.10    -1.34952    5.06728  -0.266  0.79004   
## H3.valid1:seasonSH.11    -1.37867    4.08013  -0.338  0.73550   
## H3.valid1:seasonSH.12     0.39415    4.57904   0.086  0.93142   
## H3.valid1:seasonSH.13    -3.20202    4.12508  -0.776  0.43776   
## H3.valid1:seasonSH.14    -2.80429    3.96570  -0.707  0.47961   
## H3.valid1:seasonSH.15     1.43405    4.19364   0.342  0.73244   
## H3.valid1:seasonSH.16     0.41060    3.46223   0.119  0.90562   
## H3.valid1:seasonSH.17          NA         NA      NA       NA   
## H3.valid1:anyvac1        -0.25152    0.84449  -0.298  0.76587   
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 43.85066)
## 
##     Null deviance: 58890  on 1320  degrees of freedom
## Residual deviance: 55778  on 1272  degrees of freedom
##   (85 observations deleted due to missingness)
## AIC: 8793.3
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions NO COUNTRY AND SEASON
fit3 = glm(symdur ~ ns(age, knots = 65) + anyav + anydx  + H3.valid*anyvac, data = train)
summary(fit3)
```

```
## 
## Call:
## glm(formula = symdur ~ ns(age, knots = 65) + anyav + anydx + 
##     H3.valid * anyvac, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -7.770   -2.977   -1.440    0.848  143.948  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            6.6344     0.7421   8.941  < 2e-16 ***
## ns(age, knots = 65)1   1.1954     1.2862   0.929  0.35282    
## ns(age, knots = 65)2  -0.7615     0.9899  -0.769  0.44184    
## anyav1                -1.5372     0.4952  -3.104  0.00195 ** 
## anydx1                 0.3041     0.4836   0.629  0.52956    
## H3.valid1             -0.7373     0.4788  -1.540  0.12384    
## anyvac1                0.3199     0.6026   0.531  0.59563    
## H3.valid1:anyvac1      0.6694     0.7837   0.854  0.39322    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 44.27461)
## 
##     Null deviance: 58890  on 1320  degrees of freedom
## Residual deviance: 58133  on 1313  degrees of freedom
##   (85 observations deleted due to missingness)
## AIC: 8765.9
## 
## Number of Fisher Scoring iterations: 2
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, NO COUNTRY AND SEASON
fit4 = glm(symdur ~ ns(age, knots = 65) + p.protected + anyav + anydx + H3.valid*anyvac, data = train)
summary(fit4)
```

```
## 
## Call:
## glm(formula = symdur ~ ns(age, knots = 65) + p.protected + anyav + 
##     anydx + H3.valid * anyvac, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
##  -7.647   -2.995   -1.465    0.899  144.181  
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            6.8265     0.7880   8.664  < 2e-16 ***
## ns(age, knots = 65)1   1.2408     1.2879   0.963  0.33551    
## ns(age, knots = 65)2  -0.8829     1.0041  -0.879  0.37938    
## p.protected           -0.3753     0.5170  -0.726  0.46797    
## anyav1                -1.5488     0.4956  -3.125  0.00182 ** 
## anydx1                 0.3138     0.4838   0.649  0.51677    
## H3.valid1             -0.8584     0.5071  -1.693  0.09075 .  
## anyvac1                0.3846     0.6093   0.631  0.52803    
## H3.valid1:anyvac1      0.5676     0.7963   0.713  0.47614    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 44.29056)
## 
##     Null deviance: 58890  on 1320  degrees of freedom
## Residual deviance: 58109  on 1312  degrees of freedom
##   (85 observations deleted due to missingness)
## AIC: 8767.4
## 
## Number of Fisher Scoring iterations: 2
```

```r
## Plot the fits for each model side by side
par(mfrow = c(1, 2), las = 3, mar = c(6, 6, 6, 3), cex.axis = .7)
library(gam)
plot.gam(fit00, se = T, col = 'darkslateblue', main = 'Simplified Baseline')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png)

```r
plot.gam(fit0, se = T, col = 'darkslateblue', main = 'Baseline')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-3.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-4.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-5.png)

```r
plot.gam(fit1, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-6.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-7.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-8.png)

```r
plot.gam(fit2, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-9.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-10.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-11.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-12.png)

```r
plot.gam(fit3, se = T, col = 'darkslateblue', main = 'Simple baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-13.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-14.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-15.png)

```r
plot.gam(fit4, se = T, col = 'darkslateblue', main = 'Simple baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-16.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-17.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-18.png)

```r
pdf('003Symdur_constrained.pdf')
plot.gam(fit4, se = T, col = 'darkslateblue', main = 'Simple baseline + age + imp')
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
# Predict probs for each patient
pred00 = predict(fit00, newdata = test)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-19.png)

```r
pred0 = predict(fit0, newdata = test)
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```r
pred1 = predict(fit1, newdata = test)
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```r
pred2 = predict(fit2, newdata = test)
```

```
## Warning in predict.lm(object, newdata, se.fit, scale = 1, type =
## ifelse(type == : prediction from a rank-deficient fit may be misleading
```

```r
pred3 = predict(fit3, newdata = test)
pred4 = predict(fit4, newdata = test)
  
# Estimate the simulated error rate
MSE = numeric(6); names(MSE) = c('simple.baseline', 'baseline', 'baeline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')
MSE[1] = mean((pred00 - test$symdur)^2)
MSE[2] = mean((pred0 - test$symdur)^2)
MSE[3] = mean((pred1 - test$symdur)^2)
MSE[4] = mean((pred2 - test$symdur)^2)
MSE[5] = mean((pred3 - test$symdur)^2)
MSE[6] = mean((pred4 - test$symdur)^2)
sort(MSE)
```

```
## named numeric(0)
```

```r
AIC = c(fit00$aic, fit0$aic, fit1$aic, fit2$aic, fit3$aic, fit4$aic); names(AIC) =  c('simple.baseline', 'baseline', 'baeline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')
del.AIC = sort(AIC - min(AIC))
del.AIC
```

```
##         simple.baseline     simple.baseline+age simple.baseline+age+imp 
##                0.000000                2.494667                3.964094 
##                baseline             baeline+age        baseline+age+imp 
##               26.132208               28.018619               29.877031
```

```r
load('master.AIC.table.constrained.RData')
master.AIC.table['003Symdur',  c('simple.baseline', 'baseline', 'baseline+age', 'baseline+age+imp', 'simple.baseline+age', 'simple.baseline+age+imp')] = AIC - min(AIC)
save(master.AIC.table, file = 'master.AIC.table.constrained.RData')
rm(master.AIC.table)
```




