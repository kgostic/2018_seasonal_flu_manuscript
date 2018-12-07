---
title: "002 incidence"
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


#### Visualize the relationship between age and symptoms >= 14d, given infection
![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)



#### Visualize the relationship between age and prob of H3N2 infection
![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

#### Visualize the relationship between country and prob infection

H1N1
![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

H3N2
![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)


#### Visualize the relationship between season and prob infection

H1N1
![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

H3N2
![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)




* age
* imprinting status
* vaccination status
* antiviral use
* season
* country





### Clean data for fitting

```
## 
## Argentina Australia   Austria   Belgium     Chile   Denmark   Estonia 
##      2277        79        15       820        39       141       150 
##   Germany    Greece     Japan      Peru    Poland  Portugal     Spain 
##       313       395        21       387       250         6       112 
##  Thailand        UK       USA 
##      2510        74       817
```

```
## [1] 8406
```

```
## 
## Argentina Australia   Austria   Belgium     Chile   Denmark   Estonia 
##      2516        68         6       851        25       120       108 
##   Germany    Greece     Japan      Peru    Poland  Portugal     Spain 
##       285       244        29       414       227         6        87 
##  Thailand        UK       USA 
##      2689        72       738
```

```
## [1] 8485
```






### First fit a model where the only independent variable is age
Treat vaccination, antivial use, underlying symptoms, country and season as blocking variables with fixed effects

```r
library(splines)
library(ROCR)

## Fit to country, season and medical history only (no effect of age or imprinting)
fit0 = glm(infected ~ anyav + anydx + country + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit0)
```

```
## 
## Call:
## glm(formula = infected ~ anyav + anydx + country + H3.valid * 
##     season + H3.valid * anyvac, family = binomial, data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.2466  -0.7410  -0.4540  -0.0003   3.2107  
## 
## Coefficients:
##                           Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                0.66899    0.12875   5.196 2.04e-07 ***
## anyav1                     0.88160    0.07817  11.277  < 2e-16 ***
## anydx1                    -0.22342    0.04799  -4.655 3.24e-06 ***
## countryAustralia          -0.53373    0.23125  -2.308 0.020998 *  
## countryBelgium            -0.13824    0.10795  -1.281 0.200317    
## countryChile               2.83024    0.50209   5.637 1.73e-08 ***
## countryDenmark            -0.52255    0.18417  -2.837 0.004549 ** 
## countryEstonia            -0.08143    0.17331  -0.470 0.638464    
## countryGermany            -1.10494    0.16544  -6.679 2.41e-11 ***
## countryGreece             -0.73187    0.14305  -5.116 3.12e-07 ***
## countryJapan              -0.50757    0.79388  -0.639 0.522594    
## countryOther              -0.50006    0.29630  -1.688 0.091471 .  
## countryPeru               -1.14429    0.12091  -9.464  < 2e-16 ***
## countryPoland             -0.56782    0.14902  -3.810 0.000139 ***
## countrySpain              -1.08905    0.20294  -5.366 8.04e-08 ***
## countryThailand           -0.78640    0.07709 -10.201  < 2e-16 ***
## countryUK                 -1.54112    0.31312  -4.922 8.57e-07 ***
## countryUSA                -0.71982    0.11418  -6.304 2.90e-10 ***
## H3.valid1                -16.48012  127.72199  -0.129 0.897333    
## seasonNH.10.11            -0.40795    0.18037  -2.262 0.023711 *  
## seasonNH.11.12            -5.09781    1.01008  -5.047 4.49e-07 ***
## seasonNH.12.13            -1.00311    0.15583  -6.437 1.22e-10 ***
## seasonNH.13.14            -1.26684    0.12662 -10.005  < 2e-16 ***
## seasonNH.14.15            -2.56429    0.18402 -13.935  < 2e-16 ***
## seasonNH.15.16            -0.83369    0.11981  -6.958 3.44e-12 ***
## seasonNH.16.17            -0.85965    0.34300  -2.506 0.012202 *  
## seasonSH.10               -0.19613    0.18297  -1.072 0.283769    
## seasonSH.11               -3.95916    0.33519 -11.812  < 2e-16 ***
## seasonSH.12               -2.20907    0.25602  -8.628  < 2e-16 ***
## seasonSH.13               -1.05913    0.14766  -7.173 7.36e-13 ***
## seasonSH.14               -3.23511    0.22642 -14.288  < 2e-16 ***
## seasonSH.15               -3.36797    0.21367 -15.763  < 2e-16 ***
## seasonSH.16                0.87342    0.16226   5.383 7.33e-08 ***
## anyvac1                   -0.74684    0.08859  -8.430  < 2e-16 ***
## anyvac2                    0.23809    0.55114   0.432 0.665744    
## H3.valid1:seasonNH.10.11  14.67943  127.72242   0.115 0.908499    
## H3.valid1:seasonNH.11.12  21.08773  127.72601   0.165 0.868864    
## H3.valid1:seasonNH.12.13  16.48688  127.72211   0.129 0.897291    
## H3.valid1:seasonNH.13.14  15.92597  127.72205   0.125 0.900767    
## H3.valid1:seasonNH.14.15  18.23295  127.72210   0.143 0.886484    
## H3.valid1:seasonNH.15.16  15.30726  127.72207   0.120 0.904603    
## H3.valid1:seasonNH.16.17  17.63601  127.72258   0.138 0.890177    
## H3.valid1:seasonSH.10     14.83504  127.72228   0.116 0.907533    
## H3.valid1:seasonSH.11     19.45741  127.72240   0.152 0.878918    
## H3.valid1:seasonSH.12     17.16184  127.72229   0.134 0.893111    
## H3.valid1:seasonSH.13     15.54380  127.72207   0.122 0.903136    
## H3.valid1:seasonSH.14     18.72754  127.72216   0.147 0.883426    
## H3.valid1:seasonSH.15     18.73128  127.72213   0.147 0.883403    
## H3.valid1:seasonSH.16     15.26468  127.72214   0.120 0.904868    
## H3.valid1:anyvac1          0.46474    0.11202   4.149 3.34e-05 ***
## H3.valid1:anyvac2         -0.09192    0.81304  -0.113 0.909987    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 16552  on 14890  degrees of freedom
## Residual deviance: 13598  on 14840  degrees of freedom
## AIC: 13700
## 
## Number of Fisher Scoring iterations: 15
```

```r
## H1N1 fit to age only, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit1 = glm(infected ~ anyav + anydx + country + ns(age, knots = 65) + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit1)
```

```
## 
## Call:
## glm(formula = infected ~ anyav + anydx + country + ns(age, knots = 65) + 
##     H3.valid * season + H3.valid * anyvac, family = binomial, 
##     data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.2472  -0.7416  -0.4566  -0.0003   3.2102  
## 
## Coefficients:
##                           Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                0.66961    0.13393   5.000 5.74e-07 ***
## anyav1                     0.88367    0.07821  11.298  < 2e-16 ***
## anydx1                    -0.21437    0.05002  -4.286 1.82e-05 ***
## countryAustralia          -0.53464    0.23143  -2.310 0.020881 *  
## countryBelgium            -0.13482    0.10832  -1.245 0.213248    
## countryChile               2.82660    0.50220   5.628 1.82e-08 ***
## countryDenmark            -0.52131    0.18421  -2.830 0.004655 ** 
## countryEstonia            -0.08302    0.17330  -0.479 0.631893    
## countryGermany            -1.10808    0.16562  -6.690 2.22e-11 ***
## countryGreece             -0.72889    0.14319  -5.090 3.57e-07 ***
## countryJapan              -0.51334    0.79395  -0.647 0.517914    
## countryOther              -0.50509    0.29640  -1.704 0.088361 .  
## countryPeru               -1.14752    0.12100  -9.484  < 2e-16 ***
## countryPoland             -0.56713    0.14903  -3.806 0.000141 ***
## countrySpain              -1.07570    0.20356  -5.285 1.26e-07 ***
## countryThailand           -0.78318    0.07723 -10.141  < 2e-16 ***
## countryUK                 -1.53631    0.31325  -4.904 9.37e-07 ***
## countryUSA                -0.72202    0.11434  -6.315 2.71e-10 ***
## ns(age, knots = 65)1      -0.10822    0.15478  -0.699 0.484428    
## ns(age, knots = 65)2      -0.15450    0.24103  -0.641 0.521505    
## H3.valid1                -16.47919  127.71491  -0.129 0.897333    
## seasonNH.10.11            -0.40653    0.18053  -2.252 0.024333 *  
## seasonNH.11.12            -5.09332    1.01015  -5.042 4.60e-07 ***
## seasonNH.12.13            -0.99981    0.15626  -6.399 1.57e-10 ***
## seasonNH.13.14            -1.26293    0.12714  -9.933  < 2e-16 ***
## seasonNH.14.15            -2.56139    0.18427 -13.900  < 2e-16 ***
## seasonNH.15.16            -0.82930    0.12039  -6.889 5.64e-12 ***
## seasonNH.16.17            -0.85424    0.34345  -2.487 0.012875 *  
## seasonSH.10               -0.19788    0.18301  -1.081 0.279570    
## seasonSH.11               -3.95434    0.33535 -11.792  < 2e-16 ***
## seasonSH.12               -2.20311    0.25621  -8.599  < 2e-16 ***
## seasonSH.13               -1.05487    0.14806  -7.125 1.04e-12 ***
## seasonSH.14               -3.22933    0.22678 -14.240  < 2e-16 ***
## seasonSH.15               -3.36095    0.21402 -15.704  < 2e-16 ***
## seasonSH.16                0.87875    0.16266   5.402 6.57e-08 ***
## anyvac1                   -0.74172    0.08895  -8.339  < 2e-16 ***
## anyvac2                    0.23693    0.55148   0.430 0.667475    
## H3.valid1:seasonNH.10.11  14.67873  127.71534   0.115 0.908498    
## H3.valid1:seasonNH.11.12  21.08697  127.71893   0.165 0.868862    
## H3.valid1:seasonNH.12.13  16.48712  127.71503   0.129 0.897284    
## H3.valid1:seasonNH.13.14  15.92515  127.71497   0.125 0.900767    
## H3.valid1:seasonNH.14.15  18.23403  127.71502   0.143 0.886471    
## H3.valid1:seasonNH.15.16  15.30709  127.71499   0.120 0.904599    
## H3.valid1:seasonNH.16.17  17.63639  127.71550   0.138 0.890168    
## H3.valid1:seasonSH.10     14.83621  127.71520   0.116 0.907521    
## H3.valid1:seasonSH.11     19.45704  127.71532   0.152 0.878913    
## H3.valid1:seasonSH.12     17.16027  127.71521   0.134 0.893115    
## H3.valid1:seasonSH.13     15.54372  127.71499   0.122 0.903132    
## H3.valid1:seasonSH.14     18.72631  127.71508   0.147 0.883427    
## H3.valid1:seasonSH.15     18.72966  127.71505   0.147 0.883407    
## H3.valid1:seasonSH.16     15.26402  127.71506   0.120 0.904866    
## H3.valid1:anyvac1          0.46441    0.11203   4.145 3.39e-05 ***
## H3.valid1:anyvac2         -0.08804    0.81327  -0.108 0.913793    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 16552  on 14890  degrees of freedom
## Residual deviance: 13597  on 14838  degrees of freedom
## AIC: 13703
## 
## Number of Fisher Scoring iterations: 15
```

```r
## H1N1 fit to age plus imprinting, plus fixed effects of vaccination, av use, underlying conditions, country and season
fit2 = glm(infected ~ anyav + anydx + country + ns(age, knots = 65) + p.protected + H3.valid*season + H3.valid*anyvac, data = train, family = binomial)
summary(fit2)
```

```
## 
## Call:
## glm(formula = infected ~ anyav + anydx + country + ns(age, knots = 65) + 
##     p.protected + H3.valid * season + H3.valid * anyvac, family = binomial, 
##     data = train)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.2588  -0.7429  -0.4582  -0.0003   3.1995  
## 
## Coefficients:
##                           Estimate Std. Error z value Pr(>|z|)    
## (Intercept)                0.72415    0.13710   5.282 1.28e-07 ***
## anyav1                     0.88544    0.07823  11.319  < 2e-16 ***
## anydx1                    -0.21338    0.05002  -4.266 1.99e-05 ***
## countryAustralia          -0.53490    0.23140  -2.312 0.020800 *  
## countryBelgium            -0.13795    0.10834  -1.273 0.202910    
## countryChile               2.82345    0.50290   5.614 1.97e-08 ***
## countryDenmark            -0.52435    0.18418  -2.847 0.004415 ** 
## countryEstonia            -0.09243    0.17347  -0.533 0.594162    
## countryGermany            -1.11378    0.16563  -6.725 1.76e-11 ***
## countryGreece             -0.73103    0.14324  -5.103 3.34e-07 ***
## countryJapan              -0.51105    0.79325  -0.644 0.519415    
## countryOther              -0.50892    0.29650  -1.716 0.086085 .  
## countryPeru               -1.15277    0.12104  -9.524  < 2e-16 ***
## countryPoland             -0.56953    0.14907  -3.821 0.000133 ***
## countrySpain              -1.07201    0.20329  -5.273 1.34e-07 ***
## countryThailand           -0.78745    0.07728 -10.190  < 2e-16 ***
## countryUK                 -1.53936    0.31325  -4.914 8.91e-07 ***
## countryUSA                -0.72914    0.11444  -6.371 1.87e-10 ***
## ns(age, knots = 65)1      -0.14783    0.15638  -0.945 0.344500    
## ns(age, knots = 65)2      -0.20854    0.24331  -0.857 0.391386    
## p.protected               -0.11578    0.06165  -1.878 0.060377 .  
## H3.valid1                -16.46296  127.72135  -0.129 0.897439    
## seasonNH.10.11            -0.40008    0.18064  -2.215 0.026777 *  
## seasonNH.11.12            -5.08637    1.01015  -5.035 4.77e-07 ***
## seasonNH.12.13            -0.99686    0.15635  -6.376 1.82e-10 ***
## seasonNH.13.14            -1.25598    0.12719  -9.875  < 2e-16 ***
## seasonNH.14.15            -2.56306    0.18430 -13.907  < 2e-16 ***
## seasonNH.15.16            -0.82910    0.12043  -6.885 5.79e-12 ***
## seasonNH.16.17            -0.84780    0.34380  -2.466 0.013664 *  
## seasonSH.10               -0.20264    0.18299  -1.107 0.268128    
## seasonSH.11               -3.95266    0.33530 -11.788  < 2e-16 ***
## seasonSH.12               -2.20465    0.25629  -8.602  < 2e-16 ***
## seasonSH.13               -1.05466    0.14811  -7.121 1.07e-12 ***
## seasonSH.14               -3.22560    0.22688 -14.217  < 2e-16 ***
## seasonSH.15               -3.35754    0.21403 -15.687  < 2e-16 ***
## seasonSH.16                0.87635    0.16279   5.383 7.32e-08 ***
## anyvac1                   -0.72827    0.08922  -8.162 3.29e-16 ***
## anyvac2                    0.24136    0.55217   0.437 0.662037    
## H3.valid1:seasonNH.10.11  14.67178  127.72179   0.115 0.908546    
## H3.valid1:seasonNH.11.12  21.07766  127.72537   0.165 0.868926    
## H3.valid1:seasonNH.12.13  16.48515  127.72147   0.129 0.897301    
## H3.valid1:seasonNH.13.14  15.91664  127.72141   0.125 0.900824    
## H3.valid1:seasonNH.14.15  18.23888  127.72147   0.143 0.886447    
## H3.valid1:seasonNH.15.16  15.31008  127.72143   0.120 0.904585    
## H3.valid1:seasonNH.16.17  17.62471  127.72194   0.138 0.890246    
## H3.valid1:seasonSH.10     14.84040  127.72165   0.116 0.907499    
## H3.valid1:seasonSH.11     19.45068  127.72176   0.152 0.878959    
## H3.valid1:seasonSH.12     17.16256  127.72166   0.134 0.893106    
## H3.valid1:seasonSH.13     15.54190  127.72144   0.122 0.903148    
## H3.valid1:seasonSH.14     18.71952  127.72152   0.147 0.883475    
## H3.valid1:seasonSH.15     18.72290  127.72149   0.147 0.883454    
## H3.valid1:seasonSH.16     15.26692  127.72150   0.120 0.904853    
## H3.valid1:anyvac1          0.44300    0.11262   3.934 8.37e-05 ***
## H3.valid1:anyvac2         -0.10566    0.81383  -0.130 0.896705    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for binomial family taken to be 1)
## 
##     Null deviance: 16552  on 14890  degrees of freedom
## Residual deviance: 13594  on 14837  degrees of freedom
## AIC: 13702
## 
## Number of Fisher Scoring iterations: 15
```

```r
## Plot the fits for each model side by side
par(mfrow = c(1, 2), las = 3, mar = c(6, 6, 6, 3), cex.axis = .7)
library(gam)
plot.gam(fit0, se = T, col = 'darkslateblue', main = 'Baseline')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-3.png)

```r
plot.gam(fit1, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-4.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-5.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-6.png)

```r
plot.gam(fit2, se = T, col = 'darkslateblue', main = 'Baseline + age + imp')
```

```
## Warning in preplot.gam(x, terms = terms): No terms saved for "a:b" style
## interaction terms
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-7.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-8.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-9.png)![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-10.png)

```r
pdf('002Incidence_constrained.pdf')
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

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-11.png)

```r
probs1 = predict(fit1, newdata = test, type = 'response')
probs2 = predict(fit2, newdata = test, type = 'response')

# Create a predict object
pr0 = prediction(probs0, test$infected)
pr1 = prediction(probs1, test$infected)
pr2 = prediction(probs2, test$infected)

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
pdf('002Incidence_AUC_constrained.pdf')
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
##         baseline baseline+age+imp      baeline+age 
##         0.000000         1.830317         3.362118
```

```r
load('master.AIC.table.constrained.RData')
master.AIC.table['002Incidence', c('baseline', 'baseline+age', 'baseline+age+imp')] = AIC - min(AIC)
save(master.AIC.table, file = 'master.AIC.table.constrained.RData')
rm(master.AIC.table)
```


