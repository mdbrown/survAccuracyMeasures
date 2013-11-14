survAccuracyMeasures
=============================================

This R package computes non-parametric and semi-parametric estimates of accuracy measures for risk prediction markers from survival data. It consists of the function `survAM.estimate` which estimates the *AUC*, *TPR( c )*, *FPR( c )*, *PPV( c )*, and *NPV( c )* for for a specific timepoint and marker cutoff value c. Bootstrap standard errors, and confidence intervals are also computed. 

For detailed information regarding estimation methods, see references below. 



## Tutorial



```r
# install the package from github download the package from github
if (!require("devtools")) install.packages("devtools")
devtools::install_github("survMarkerTwoPhase", "mdbrown")
```



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



Estimate all measures, using the bootstrap to estimate standard errors. First we obtain non-parametric estimates using double inverse probablity weighting:


```r
# Estimate all the measures at future time 2, with a marker cutpoint at 0.

# non-parametric estimates
survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    ESTmethod = "NP", predict.time = 1, cutpoint = 0, bootstraps = 50)
```

```
## 
## Non-Parametric estimates of accuracy measures:
##    (SE's calculated using the bootstrap)
## 
##         estimate     se      lower 0.95  upper 0.95
## AUC        0.786     0.026         0.730       0.833 
## TPR(c)     0.769     0.047         0.665       0.848 
## FPR(c)     0.410     0.029         0.356       0.467 
## PPV(c)     0.231     0.035         0.169       0.306 
## NPV(c)     0.941     0.014         0.908       0.963 
## 
##  marker cutpoint: c = 0
```


Alternatively, we can calculate semi-parametric estimates based on a proportional hazards model:



```r
# semi-parametric estimates assuming a cox model
survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    ESTmethod = "SP", predict.time = 1, cutpoint = 0, bootstraps = 50)
```

```
## 
## Semi-Parametric estimates of accuracy measures:
##    (SE's calculated using the bootstrap)
## 
##         estimate     se      lower 0.95  upper 0.95
## coef       1.010     0.077         0.859       1.161 
## AUC        0.768     0.020         0.728       0.805 
## TPR(c)     0.788     0.026         0.732       0.835 
## FPR(c)     0.412     0.028         0.359       0.468 
## PPV(c)     0.241     0.025         0.196       0.294 
## NPV(c)     0.943     0.007         0.928       0.955 
## 
##  marker cutpoint: c = 0
```




Only estimate the $AUC$ and $TPR(0)$. 


```r
tmp <- survAM.estimate(time = SimData$survTime, event = SimData$status, marker = SimData$Y, 
    predict.time = 2, measures = c("AUC", "TPR"), bootstraps = 50, cutpoint = 0)
tmp
```

```
## 
## Non-Parametric estimates of accuracy measures:
##    (SE's calculated using the bootstrap)
## 
##         estimate     se      lower 0.95  upper 0.95
## AUC        0.750     0.027         0.693       0.800 
## TPR(c)     0.699     0.050         0.593       0.787 
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
## upper 0.7996 0.7874
## lower 0.6932 0.5929
```


For more information see `?survAM.estimate`. 


### References
Liu D, Cai T, Zheng Y. Evaluating the predictive value of biomarkers with stratified case-cohort design. *Biometrics* 2012, 4: 1219-1227.

Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance of markers for future events. *Lifetime Data Analysis.* 2008, 14: 86-113.

Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent predictive values of prognostic biomarkers with failure time outcome. *JASA* 2008, 103: 362-368.














