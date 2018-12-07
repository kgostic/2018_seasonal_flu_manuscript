---
title: "002 prolonged symptoms"
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


Notes on duration of symptoms:
A review by Carrat et al. (2007) of challenge studies in healthy human volunteers suggests that influenza symptoms typically resolve fully within 8-9 days of inoculation. This led us to classify patients with the following characteristics as having long-lasting symptoms:

* Symptoms lasting 10 or more days prior to study enrollment (this implies that symptoms on day 10 were still serious enough to send patients to the doctor, and most likely lasted at least a few days longer.)
* Symptoms not resolved on day 14 of enrollment (implies that symptoms lasted longer than 14 days, as most patients had symptoms for at least one day prior to enrollment.)


```r
# use the above definition to define a column for prolonged symptoms
dat.002$prolonged.symp = as.numeric(dat.002$symdur >= 10 | dat.002$resolved == "0")
```


#### Visualize the relationship between age and prob of H1N1 infection
![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)



#### Visualize the relationship between age and prob of H3N2 infection
![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

#### Visualize the relationship between country and prolonged symptoms

H1N1

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

H3N2

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)


#### Visualize the relationship between season and prob infection

H1N1

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

H3N2

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)



### Clean data for fitting

```
## 
## Argentina Australia   Austria   Belgium     Chile   Denmark   Estonia 
##      3368        97        15      1300        38       200       217 
##   Germany    Greece     Japan      Peru    Poland  Portugal     Spain 
##       352       397        32       522       355         6       135 
##  Thailand        UK       USA 
##      3276        79       922
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
## [1] 11311
```

```
## 
## Argentina Australia   Austria   Belgium     Chile   Denmark   Estonia 
##      3368        97        15      1300        38       200       217 
##   Germany    Greece     Japan      Peru    Poland  Portugal     Spain 
##       352       397        32       522       355         6       135 
##  Thailand        UK       USA 
##      3276        79       922
```

```
## [1] 11311
```






### First fit a model where the only independent variable is age
Treat vaccination, antivial use, underlying symptoms, country and season as blocking variables with fixed effects

```r
library(splines)
library(ROCR)

## Fit to country, season and medical history only (no effect of age or imprinting)
fit0 = glm(prolonged.symp ~ anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit0)
```

```
## 
## Call:
## glm(formula = prolonged.symp ~ anyav + anydx + country + H3.valid * 
##     season + H3.valid * anyvac, family = binomial, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.5374  -0.5319  -0.3786  -0.2852   2.8442  
## 
## Coefficients:
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)              -2.449e+00  1.420e-01 -17.249  < 2e-16 ***
## anyav1                    9.864e-02  7.552e-02   1.306 0.191497    
## anydx1                    3.040e-01  4.584e-02   6.633 3.28e-11 ***
## countryAustralia          1.570e+00  1.848e-01   8.493  < 2e-16 ***
## countryBelgium            1.158e+00  1.131e-01  10.240  < 2e-16 ***
## countryChile              1.201e+00  2.622e-01   4.579 4.68e-06 ***
## countryDenmark            1.525e+00  1.496e-01  10.191  < 2e-16 ***
## countryEstonia            1.130e-01  1.849e-01   0.611 0.541281    
## countryGermany            1.839e+00  1.280e-01  14.366  < 2e-16 ***
## countryGreece             1.302e-01  1.604e-01   0.812 0.416978    
## countryJapan              1.024e+00  2.935e-01   3.487 0.000488 ***
## countryOther             -1.220e+01  1.377e+02  -0.089 0.929375    
## countryPeru              -7.224e-01  1.443e-01  -5.008 5.50e-07 ***
## countryPoland             7.645e-01  1.390e-01   5.501 3.78e-08 ***
## countrySpain              1.242e+00  1.720e-01   7.219 5.23e-13 ***
## countryThailand          -8.315e-01  9.468e-02  -8.782  < 2e-16 ***
## countryUK                 2.709e+00  1.965e-01  13.781  < 2e-16 ***
## countryUSA                1.562e+00  1.123e-01  13.905  < 2e-16 ***
## H3.valid1                 5.575e-03  1.399e-01   0.040 0.968204    
## seasonNH.10.11            9.769e-02  1.876e-01   0.521 0.602623    
## seasonNH.11.12            4.781e-01  1.674e-01   2.855 0.004304 ** 
## seasonNH.12.13            4.328e-01  1.340e-01   3.230 0.001238 ** 
## seasonNH.13.14            4.357e-01  1.245e-01   3.500 0.000465 ***
## seasonNH.14.15            2.521e-01  1.249e-01   2.018 0.043542 *  
## seasonNH.15.16            1.120e-02  1.285e-01   0.087 0.930555    
## seasonNH.16.17           -7.250e-02  2.830e-01  -0.256 0.797843    
## seasonSH.10               3.397e-01  2.540e-01   1.337 0.181063    
## seasonSH.11               3.983e-01  1.836e-01   2.169 0.030101 *  
## seasonSH.12               5.650e-01  2.072e-01   2.727 0.006385 ** 
## seasonSH.13               4.303e-01  1.717e-01   2.506 0.012224 *  
## seasonSH.14               2.606e-01  1.608e-01   1.621 0.105049    
## seasonSH.15              -2.265e-01  1.695e-01  -1.336 0.181634    
## seasonSH.16              -7.463e-01  2.211e-01  -3.375 0.000739 ***
## seasonSH.17              -1.270e+01  8.827e+02  -0.014 0.988521    
## anyvac1                   5.106e-02  6.735e-02   0.758 0.448394    
## anyvac2                  -5.917e-01  5.296e-01  -1.117 0.263905    
## H3.valid1:seasonNH.10.11  6.246e-02  2.666e-01   0.234 0.814794    
## H3.valid1:seasonNH.11.12  3.021e-03  2.340e-01   0.013 0.989700    
## H3.valid1:seasonNH.12.13 -1.086e-01  1.905e-01  -0.570 0.568767    
## H3.valid1:seasonNH.13.14 -6.164e-02  1.754e-01  -0.351 0.725233    
## H3.valid1:seasonNH.14.15 -1.371e-02  1.769e-01  -0.078 0.938213    
## H3.valid1:seasonNH.15.16 -1.699e-02  1.820e-01  -0.093 0.925611    
## H3.valid1:seasonNH.16.17  4.559e-02  4.056e-01   0.112 0.910505    
## H3.valid1:seasonSH.10    -1.784e-01  3.432e-01  -0.520 0.603079    
## H3.valid1:seasonSH.11     1.926e-02  2.182e-01   0.088 0.929668    
## H3.valid1:seasonSH.12    -7.356e-02  2.704e-01  -0.272 0.785594    
## H3.valid1:seasonSH.13    -2.043e-02  2.126e-01  -0.096 0.923419    
## H3.valid1:seasonSH.14    -2.494e-02  1.960e-01  -0.127 0.898766    
## H3.valid1:seasonSH.15     1.590e-02  2.111e-01   0.075 0.939972    
## H3.valid1:seasonSH.16     7.941e-02  2.944e-01   0.270 0.787368    
## H3.valid1:seasonSH.17    -9.504e-03  1.248e+03   0.000 0.999994    
## H3.valid1:anyvac1         3.929e-03  9.440e-02   0.042 0.966802    
## H3.valid1:anyvac2         2.310e-01  7.554e-01   0.306 0.759713    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 18567  on 21537  degrees of freedom
## Residual deviance: 15916  on 21485  degrees of freedom
## AIC: 16022
## 
## Number of Fisher Scoring iterations: 13
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(prolonged.symp ~ anyav + anydx + country + ns(age, knots = 65) + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit1)
```

```
## 
## Call:
## glm(formula = prolonged.symp ~ anyav + anydx + country + ns(age, 
##     knots = 65) + H3.valid * season + H3.valid * anyvac, family = binomial, 
##     data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7311  -0.5462  -0.3831  -0.2789   2.9425  
## 
## Coefficients:
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)              -2.638e+00  1.472e-01 -17.915  < 2e-16 ***
## anyav1                    8.951e-02  7.569e-02   1.183 0.236940    
## anydx1                    1.963e-01  4.773e-02   4.113 3.91e-05 ***
## countryAustralia          1.563e+00  1.862e-01   8.392  < 2e-16 ***
## countryBelgium            1.116e+00  1.134e-01   9.844  < 2e-16 ***
## countryChile              1.234e+00  2.633e-01   4.685 2.79e-06 ***
## countryDenmark            1.518e+00  1.499e-01  10.127  < 2e-16 ***
## countryEstonia            1.339e-01  1.850e-01   0.723 0.469391    
## countryGermany            1.867e+00  1.285e-01  14.533  < 2e-16 ***
## countryGreece             8.873e-02  1.605e-01   0.553 0.580461    
## countryJapan              1.095e+00  2.943e-01   3.720 0.000200 ***
## countryOther             -1.213e+01  1.377e+02  -0.088 0.929797    
## countryPeru              -7.088e-01  1.445e-01  -4.906 9.29e-07 ***
## countryPoland             7.644e-01  1.391e-01   5.495 3.91e-08 ***
## countrySpain              1.164e+00  1.723e-01   6.758 1.40e-11 ***
## countryThailand          -8.739e-01  9.499e-02  -9.200  < 2e-16 ***
## countryUK                 2.673e+00  1.973e-01  13.545  < 2e-16 ***
## countryUSA                1.621e+00  1.127e-01  14.378  < 2e-16 ***
## ns(age, knots = 65)1      1.249e+00  1.478e-01   8.451  < 2e-16 ***
## ns(age, knots = 65)2      7.294e-01  2.154e-01   3.387 0.000707 ***
## H3.valid1                 4.145e-03  1.402e-01   0.030 0.976418    
## seasonNH.10.11            4.319e-02  1.883e-01   0.229 0.818588    
## seasonNH.11.12            3.930e-01  1.684e-01   2.334 0.019606 *  
## seasonNH.12.13            3.422e-01  1.351e-01   2.533 0.011315 *  
## seasonNH.13.14            3.699e-01  1.252e-01   2.955 0.003126 ** 
## seasonNH.14.15            1.722e-01  1.257e-01   1.371 0.170436    
## seasonNH.15.16           -5.678e-02  1.291e-01  -0.440 0.659984    
## seasonNH.16.17           -1.981e-01  2.848e-01  -0.695 0.486845    
## seasonSH.10               3.463e-01  2.560e-01   1.353 0.176098    
## seasonSH.11               3.296e-01  1.841e-01   1.790 0.073448 .  
## seasonSH.12               5.175e-01  2.077e-01   2.492 0.012714 *  
## seasonSH.13               3.652e-01  1.723e-01   2.119 0.034079 *  
## seasonSH.14               1.838e-01  1.613e-01   1.140 0.254457    
## seasonSH.15              -3.083e-01  1.701e-01  -1.813 0.069831 .  
## seasonSH.16              -8.298e-01  2.218e-01  -3.742 0.000183 ***
## seasonSH.17              -1.244e+01  8.827e+02  -0.014 0.988752    
## anyvac1                  -2.492e-02  6.821e-02  -0.365 0.714887    
## anyvac2                  -6.004e-01  5.275e-01  -1.138 0.255056    
## H3.valid1:seasonNH.10.11  6.479e-02  2.673e-01   0.242 0.808490    
## H3.valid1:seasonNH.11.12  4.728e-03  2.347e-01   0.020 0.983926    
## H3.valid1:seasonNH.12.13 -1.079e-01  1.912e-01  -0.564 0.572551    
## H3.valid1:seasonNH.13.14 -6.195e-02  1.757e-01  -0.353 0.724371    
## H3.valid1:seasonNH.14.15 -1.120e-02  1.773e-01  -0.063 0.949628    
## H3.valid1:seasonNH.15.16 -1.566e-02  1.823e-01  -0.086 0.931565    
## H3.valid1:seasonNH.16.17  4.483e-02  4.076e-01   0.110 0.912433    
## H3.valid1:seasonSH.10    -1.798e-01  3.453e-01  -0.521 0.602580    
## H3.valid1:seasonSH.11     1.836e-02  2.187e-01   0.084 0.933112    
## H3.valid1:seasonSH.12    -7.053e-02  2.710e-01  -0.260 0.794629    
## H3.valid1:seasonSH.13    -2.290e-02  2.129e-01  -0.108 0.914356    
## H3.valid1:seasonSH.14    -2.475e-02  1.964e-01  -0.126 0.899733    
## H3.valid1:seasonSH.15     1.670e-02  2.114e-01   0.079 0.937033    
## H3.valid1:seasonSH.16     8.072e-02  2.949e-01   0.274 0.784335    
## H3.valid1:seasonSH.17    -8.436e-03  1.248e+03   0.000 0.999995    
## H3.valid1:anyvac1         4.292e-03  9.477e-02   0.045 0.963878    
## H3.valid1:anyvac2         2.231e-01  7.524e-01   0.297 0.766788    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 18567  on 21537  degrees of freedom
## Residual deviance: 15845  on 21483  degrees of freedom
## AIC: 15955
## 
## Number of Fisher Scoring iterations: 13
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(prolonged.symp ~ anyav + anydx + country + ns(age, knots = 65) + p.protected + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit2)
```

```
## 
## Call:
## glm(formula = prolonged.symp ~ anyav + anydx + country + ns(age, 
##     knots = 65) + p.protected + H3.valid * season + H3.valid * 
##     anyvac, family = binomial, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -1.7300  -0.5467  -0.3830  -0.2789   2.9423  
## 
## Coefficients:
##                            Estimate Std. Error z value Pr(>|z|)    
## (Intercept)              -2.635e+00  1.491e-01 -17.668  < 2e-16 ***
## anyav1                    8.950e-02  7.569e-02   1.183 0.236992    
## anydx1                    1.963e-01  4.773e-02   4.113 3.91e-05 ***
## countryAustralia          1.563e+00  1.862e-01   8.392  < 2e-16 ***
## countryBelgium            1.116e+00  1.134e-01   9.844  < 2e-16 ***
## countryChile              1.234e+00  2.634e-01   4.685 2.80e-06 ***
## countryDenmark            1.518e+00  1.499e-01  10.128  < 2e-16 ***
## countryEstonia            1.339e-01  1.850e-01   0.724 0.469371    
## countryGermany            1.867e+00  1.285e-01  14.533  < 2e-16 ***
## countryGreece             8.872e-02  1.605e-01   0.553 0.580504    
## countryJapan              1.095e+00  2.943e-01   3.719 0.000200 ***
## countryOther             -1.213e+01  1.377e+02  -0.088 0.929797    
## countryPeru              -7.088e-01  1.445e-01  -4.906 9.28e-07 ***
## countryPoland             7.644e-01  1.391e-01   5.495 3.91e-08 ***
## countrySpain              1.164e+00  1.723e-01   6.758 1.39e-11 ***
## countryThailand          -8.739e-01  9.499e-02  -9.201  < 2e-16 ***
## countryUK                 2.673e+00  1.973e-01  13.544  < 2e-16 ***
## countryUSA                1.621e+00  1.127e-01  14.378  < 2e-16 ***
## ns(age, knots = 65)1      1.250e+00  1.479e-01   8.451  < 2e-16 ***
## ns(age, knots = 65)2      7.297e-01  2.154e-01   3.388 0.000704 ***
## p.protected              -6.941e-03  5.595e-02  -0.124 0.901268    
## H3.valid1                 4.985e-03  1.404e-01   0.036 0.971672    
## seasonNH.10.11            4.345e-02  1.883e-01   0.231 0.817556    
## seasonNH.11.12            3.933e-01  1.684e-01   2.335 0.019519 *  
## seasonNH.12.13            3.424e-01  1.351e-01   2.534 0.011267 *  
## seasonNH.13.14            3.701e-01  1.252e-01   2.956 0.003113 ** 
## seasonNH.14.15            1.723e-01  1.257e-01   1.371 0.170398    
## seasonNH.15.16           -5.686e-02  1.291e-01  -0.441 0.659489    
## seasonNH.16.17           -1.975e-01  2.848e-01  -0.693 0.488043    
## seasonSH.10               3.465e-01  2.560e-01   1.354 0.175869    
## seasonSH.11               3.300e-01  1.841e-01   1.792 0.073134 .  
## seasonSH.12               5.176e-01  2.077e-01   2.492 0.012693 *  
## seasonSH.13               3.654e-01  1.723e-01   2.120 0.033968 *  
## seasonSH.14               1.841e-01  1.613e-01   1.142 0.253655    
## seasonSH.15              -3.080e-01  1.701e-01  -1.811 0.070209 .  
## seasonSH.16              -8.296e-01  2.218e-01  -3.741 0.000184 ***
## seasonSH.17              -1.244e+01  8.827e+02  -0.014 0.988752    
## anyvac1                  -2.408e-02  6.853e-02  -0.351 0.725278    
## anyvac2                  -6.002e-01  5.275e-01  -1.138 0.255212    
## H3.valid1:seasonNH.10.11  6.424e-02  2.674e-01   0.240 0.810126    
## H3.valid1:seasonNH.11.12  4.095e-03  2.347e-01   0.017 0.986082    
## H3.valid1:seasonNH.12.13 -1.084e-01  1.912e-01  -0.567 0.570812    
## H3.valid1:seasonNH.13.14 -6.234e-02  1.757e-01  -0.355 0.722750    
## H3.valid1:seasonNH.14.15 -1.128e-02  1.773e-01  -0.064 0.949269    
## H3.valid1:seasonNH.15.16 -1.552e-02  1.823e-01  -0.085 0.932179    
## H3.valid1:seasonNH.16.17  4.358e-02  4.078e-01   0.107 0.914884    
## H3.valid1:seasonSH.10    -1.802e-01  3.454e-01  -0.522 0.601786    
## H3.valid1:seasonSH.11     1.750e-02  2.188e-01   0.080 0.936256    
## H3.valid1:seasonSH.12    -7.079e-02  2.710e-01  -0.261 0.793891    
## H3.valid1:seasonSH.13    -2.341e-02  2.130e-01  -0.110 0.912461    
## H3.valid1:seasonSH.14    -2.550e-02  1.965e-01  -0.130 0.896761    
## H3.valid1:seasonSH.15     1.587e-02  2.115e-01   0.075 0.940180    
## H3.valid1:seasonSH.16     8.020e-02  2.950e-01   0.272 0.785723    
## H3.valid1:seasonSH.17    -8.326e-03  1.248e+03   0.000 0.999995    
## H3.valid1:anyvac1         2.562e-03  9.579e-02   0.027 0.978658    
## H3.valid1:anyvac2         2.229e-01  7.524e-01   0.296 0.767045    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 18567  on 21537  degrees of freedom
## Residual deviance: 15845  on 21482  degrees of freedom
## AIC: 15957
## 
## Number of Fisher Scoring iterations: 13
```

```r
# Calculate 


## Plot the fits for each model side by side
par(mfrow = c(1, 2), las = 3, mar = c(6, 6, 6, 3), cex.axis = .7)
library(gam)
plot.gam(fit0, se = T, col = 'darkslateblue', main = 'Baseline')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-2.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-3.png)

```r
plot.gam(fit1, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-4.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-5.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-6.png)

```r
plot.gam(fit2, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-7.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-8.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-9.png)![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-10.png)

```r
pdf('002Duration_constrained.pdf')
plot.gam(fit2, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
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

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-11.png)

```r
probs1 = predict(fit1, newdata = test, type = 'response')
probs2 = predict(fit2, newdata = test, 'response')

# Create a predict object
pr0 = prediction(probs0, test$prolonged.symp)
pr1 = prediction(probs1, test$prolonged.symp)
pr2 = prediction(probs2, test$prolonged.symp)

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

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-12.png)

```r
pdf('002Duration_AUC_constrained.pdf')
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
# Calculate AUC for test data


AIC = c(fit0$aic, fit1$aic, fit2$aic); names(AIC) =  c('baseline', 'baeline+age', 'baseline+age+imp')
del.AIC = sort(AIC - min(AIC))
del.AIC
```

```
##      baeline+age baseline+age+imp         baseline 
##          0.00000          1.98461         66.89315
```

```r
load('master.AIC.table.constrained.RData')
master.AIC.table['002Duration', c('baseline', 'baseline+age', 'baseline+age+imp')] = AIC - min(AIC)
```

```
## Error in `[<-`(`*tmp*`, "002Duration", c("baseline", "baseline+age", "baseline+age+imp"), : subscript out of bounds
```

```r
save(master.AIC.table, file = 'master.AIC.table.constrained.RData')
rm(master.AIC.table)
```


