survAccuracyMeasures
=============================================

This R package computes semi-parametric estimates of accuracy measures for risk prediction markers from survival data. It consists of the function `survAM.estimate` which estimates the $AUC$, $TPR(c)$, $FPR(c)$, $PPV(c)$, and $NPV(c)$ for for a specific timepoint and marker cutoff value c. Standard errors, and confidence intervals are also computed. Either analytic or bootstrap standard errors can be computed. Estimation of accuracy measures for the case-cohort design is also provided.

For more information, see references below. 


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



Estimate all measures, using asymptotic normality to estimate standard errors. 


```r
# Estimate all the measures at future time 2, with a marker cutpoint at 0.

survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    predict.time = 2, cutpoint = 0, SEmethod = "normal")
```

```
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.010     0.085         0.842       1.177 
## AUC        0.775     0.019         0.736       0.809 
## TPR(c)     0.768     0.027         0.711       0.816 
## FPR(c)     0.377     0.024         0.332       0.425 
## PPV(c)     0.375     0.031         0.317       0.438 
## NPV(c)     0.901     0.013         0.872       0.924 
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
##         estimate     se      lower 0.95  upper 0.95
## coef       1.010     0.081         0.852       1.168 
## AUC        0.775     0.017         0.740       0.806 
## TPR(c)     0.768     0.027         0.711       0.817 
## 
##  marker cutpoint: c = 0
```



```r
# access the estimates
tmp$estimates
```

```
##   coef    AUC    TPR 
## 1.0099 0.7748 0.7679
```

```r

# and the confidence bounds
tmp$CIbounds
```

```
##        coef    AUC    TPR
## upper 1.168 0.8062 0.8168
## lower 0.852 0.7398 0.7105
```


For more information see `?survAM.estimate`. 



### References
Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














