survAccuracyMeasures
=============================================

This R package estimates accuracy measures for risk prediction markers from survival data. It consists of the function `survEstMeasures` which estimates the $AUC$, $TPR(c)$, $FPR(c)$, $PPV(c)$, and $NPV(c)$ for for a specific timepoint and marker cutoff value c. Standard errors, and confidence intervals are also computed. Either analytic or bootstrap standard errors can be computed.

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
## 1    2.564      0 -0.48210
## 2    1.856      0  0.42668
## 3   20.171      1 -0.79753
## 4   20.950      1  0.18762
## 5   23.277      1  0.08946
## 6    6.052      1 -0.10689
```



Estimate all measures, using asymptotic normality to estimate standard errors. 


```r
# Estimate all the measures at future time 2, with a marker cutpoint at 0.

survEstMeasures(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    predict.time = 2, cutpoint = 0, SEmethod = "normal")
```

```
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.131     0.066         1.002       1.260 
## AUC        0.793     0.015         0.763       0.820 
## TPR(c)     0.806     0.021         0.761       0.844 
## FPR(c)     0.396     0.024         0.351       0.443 
## PPV(c)     0.358     0.027         0.306       0.413 
## NPV(c)     0.919     0.010         0.898       0.937 
## 
##  marker cutpoint: c = 0
```


Only estimate the $AUC$ and $TPR(0)$. This time use bootstrapping to obtain the standard errors. 


```r
tmp <- survEstMeasures(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    predict.time = 2, measures = c("AUC", "TPR"), SEmethod = "bootstrap", bootstraps = 50, 
    cutpoint = 0)
tmp
```

```
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.131     0.070         0.993       1.268 
## AUC        0.793     0.016         0.760       0.823 
## TPR(c)     0.806     0.024         0.756       0.848 
## 
##  marker cutpoint: c = 0
```



```r
# access the estimates
tmp$estimates
```

```
##   coef    AUC    TPR 
## 1.1308 0.7934 0.8062
```

```r

# and the confidence bounds
tmp$CIbounds
```

```
##         coef    AUC    TPR
## upper 1.2683 0.8234 0.8484
## lower 0.9934 0.7597 0.7556
```


For more information see `?survEstMeasures`. 


### References
Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














