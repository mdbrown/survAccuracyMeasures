survAccuracyMeasures
=============================================

This R package estimates accuracy measures for risk prediction markers from survival data. It consists of the function `survEstMeasures` which estimates the $AUC$, $TPR(c)$, $FPR(c)$, $PPV(c)$, and $NPV(c)$ for for a specific timepoint and marker cutoff value c. Standard errors, and confidence intervals are also computed. Either analytic or bootstrap standard errors can be computed. Estimation of accuracy measures under the case-cohort design is also provided. 

For more information, see references below. 


## Tutorial



```r
library(survAccuracyMeasures)
```

```
## Loading required package: survival
```

```
## Loading required package: splines
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

survEstMeasures(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
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
tmp <- survEstMeasures(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    predict.time = 2, measures = c("AUC", "TPR"), SEmethod = "bootstrap", bootstraps = 50, 
    cutpoint = 0)
tmp
```

```
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.010     0.084         0.845       1.175 
## AUC        0.775     0.017         0.739       0.807 
## TPR(c)     0.768     0.025         0.715       0.814 
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
##         coef    AUC    TPR
## upper 1.1751 0.8068 0.8137
## lower 0.8446 0.7391 0.7148
```


For more information see `?survEstMeasures`. 

### Case-Cohort Design

Estimation using a case-cohort subcohort design is also permitted. Sample weights must first be calculated. 


```r
# generate a sub-cohort from SimData
set.seed(12321)
# create a sample index. 1 if sampled, 0 if not
N <- nrow(SimData)
sampleInd <- rep(0, N)

# sample all with observed failure time. (200 individuals)
sampleInd[SimData$status == 1] <- 1

# sample 150 more observations from the entire data set without
# replacement
sampleInd[sample(1:N, 150)] <- 1

table(sampleInd)  #total number of subcohort is 293 
```

```
## sampleInd
##   0   1 
## 207 293
```

```r

## now calculate sample weights first calculate the Pr(Sampled from
## cohort) for each observation
sampleProb <- numeric(500)
# all non-censored observations were sampled, so their sample probability
# is 1
sampleProb[SimData$status == 1] <- 1
sampleProb[SimData$status == 0] <- 150/N

SimData$weights <- 1/sampleProb

subCohortData <- SimData[sampleInd == 1, ]

# estimate accuracy measures using only the subcohort data
survEstMeasures(time = subCohortData$survTime, event = subCohortData$status, 
    marker = subCohortData$Y, weights = subCohortData$weights, predict.time = 2, 
    cutpoint = 0, SEmethod = "normal")
```

```
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.088     0.167         0.761       1.414 
## AUC        0.788     0.034         0.713       0.847 
## TPR(c)     0.770     0.058         0.639       0.864 
## FPR(c)     0.347     0.057         0.245       0.465 
## PPV(c)     0.386     0.042         0.308       0.471 
## NPV(c)     0.909     0.020         0.862       0.942 
## 
##  marker cutpoint: c = 0
```




### References
Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














