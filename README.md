survAccuracyMeasures
=============================================

This R package computes non-parametric and semi-parametric estimates of accuracy measures for risk prediction markers from survival data. It consists of the function `survAM.estimate` which estimates the *AUC*, *TPR( c )*, *FPR( c )*, *PPV( c )*, and *NPV( c )* for for a specific timepoint and marker cutoff value c. Standard errors, and confidence intervals are also computed. Either analytic or bootstrap standard errors can be computed. 

For detailed information regarding estimation methods, see references below. 



## Tutorial



```r
library(survAccuracyMeasures)
```

```
## Loading required package: survival Loading required package: splines
```

```r

# simulated data for illustration
data(SimData)
head(SimData)
```

```
##   survTime status        Y
## 1   0.1197      1  1.49310
## 2   1.0231      0 -0.73260
## 3   0.8282      0 -0.50211
## 4   2.0875      1  0.65758
## 5   4.6827      1  1.57806
## 6   0.3001      1  0.02419
```



Estimate all measures, using asymptotic normality to estimate standard errors. First we obtain non-parametric estimates using double inverse probablity weighting:


```r
# Estimate all the measures at future time 2, with a marker cutpoint at 0.

# non-parametric estimates
survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    ESTmethod = "NP", predict.time = 1, cutpoint = 0, SEmethod = "normal")
```

```
## 
## Non-Parametric estimates of accuracy measures:
##    (SE's calculated using normal approximation)
## 
##         estimate     se      lower 0.95  upper 0.95
## AUC        0.786     0.031         0.718       0.841 
## TPR(c)     0.769     0.052         0.652       0.856 
## FPR(c)     0.410     0.025         0.362       0.461 
## PPV(c)     0.231     0.030         0.178       0.294 
## NPV(c)     0.941     0.015         0.904       0.964 
## 
##  marker cutpoint: c = 0
```


Alternatively, we can calculate semi-parametric estimates based on a proportional hazards model:



```r
# semi-parametric estimates assuming a cox model
survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    ESTmethod = "SP", predict.time = 1, cutpoint = 0, SEmethod = "normal")
```

```
## 
## Semi-Parametric estimates of accuracy measures:
##    (SE's calculated using normal approximation)
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.010     0.085         0.842       1.177 
## AUC        0.768     0.020         0.727       0.805 
## TPR(c)     0.788     0.027         0.731       0.836 
## FPR(c)     0.412     0.023         0.368       0.459 
## PPV(c)     0.241     0.028         0.191       0.300 
## NPV(c)     0.943     0.009         0.924       0.958 
## 
##  marker cutpoint: c = 0
```




Only estimate the $AUC$ and $TPR(0)$. This time use bootstrapping to obtain the standard errors. 


```r
tmp <- survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    predict.time = 2, measures = c("AUC", "TPR"), SEmethod = "bootstrap", bootstraps = 50, 
    cutpoint = 0)
tmp
```

```
## 
## Non-Parametric estimates of accuracy measures:
##    (SE's calculated using the bootstrap)
## 
##         estimate     se      lower 0.95  upper 0.95
## AUC        0.750     0.026         0.696       0.798 
## TPR(c)     0.699     0.041         0.614       0.773 
## 
##  marker cutpoint: c = 0
```



```r
# access the estimates
tmp$estimates
```

```
##        AUC   TPR
## 266 0.7501 0.699
```

```r

# and the confidence bounds
tmp$CIbounds
```

```
##          AUC    TPR
## upper 0.7976 0.7726
## lower 0.6958 0.6135
```


For more information see `?survAM.estimate`. 

### Validation of R package
To validate the accuracy of estimates produced by the package, we ran several simulations under many different scenarios. Results for a single example for semi-parametric estimates (nonparametric estimates will be added soon), where *AUC = 0.75* is shown below. 1,000 cohort data sets were simulated from a proportional hazards model with sample size *n = 1,000* and *70%* censoring. Mean estimates for summary measures and SE are shown below, which can be compared to the true measure values and the empirical SE, respectively.



        | True Value | Mean(Est.) | Emp. SE | Mean(Est. SE)
--------|------------|----------|---------|-----------    
      ß |   0.879 |   0.880 | 0.064 |     0.064
    AUC |   0.750 |   0.749 | 0.015 |     0.015
 TPR(0) |   0.770 |   0.769 | 0.020 |     0.020
 FPR(0) |   0.420 |   0.420 | 0.016 |     0.017
 PPV(0) |   0.350 |   0.350 | 0.021 |     0.021
 NPV(0) |   0.900 |   0.895 | 0.010 |     0.011


### References
Liu D, Cai T, Zheng Y. Evaluating the predictive value of biomarkers with stratified case-cohort design. *Biometrics* 2012, 4: 1219-1227.

Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














