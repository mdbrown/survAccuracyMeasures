\name{survAM.estimate}
\alias{survAM.estimate}
\title{ Nonparametric and Semiparametric estimates of accuracy measures for a risk prediction marker from survival data}
\usage{
  survAM.estimate(time, event, marker,
                             data, 
                             predict.time,  
                             marker.cutpoint = 'median', 
                             estimation.method = "IPW", 
                             ci.method = "logit.transformed",
                             se.method = "bootstrap",
                             bootstraps = 1000, 
                             alpha=0.05)
}
\arguments{
  \item{time}{time to event variable}

  \item{event}{indicator for the status of
  event of interest. event = 0 for censored observations,
  and event = 1 for event of interest.}

  \item{marker}{marker variable of
  interest}
  
  \item{data}{data frame in which to look for input variables. }
  
  \item{predict.time}{numeric value of the timepoint of
  interest for which to estimate the risk measures}

  \item{marker.cutpoint}{numeric value indicating the value of the
  cutpoint 'c' at which to estimate 'FPR', 'TPR', 'NPV' or
  'PPV'. default is \code{'median'} which takes cutpoint as the marker median.}
  
  \item{estimation.method}{ Either "IPW" for non-parametric IPW estimates (default) or 
  "Cox" for semi-parametric estimates that use a Cox proportional hazards model.

  }

  \item{ci.method}{character string of either
  'logit.transformed' (default) or 'standard' indicating
  whether normal approximated confidence intervals should
  be calculated using logistic transformed values or the
  standard method.}

  \item{se.method}{ Method to calculate standard errors for estimates. Options are "bootstrap" (default) or "asymptotic". Asymptotic estimates are based on large sample calculations and will not hold in small samples. Please see referenced papers for more information.}

  \item{bootstraps}{if se.method = 'bootstrap', number of
  bootstrap replicates to use to estimate the SE.}

  \item{alpha}{alpha value for confidence intervals.
  (1-alpha)*100% confidence intervals are provided. default
  is alpha = 0.05.}
}
\value{
  a list with components \item{estimates}{point estimates
  for risk measures} \item{se}{standard errors for
  estimates} \item{CIbounds}{bounds for (1-alpha)*100%
  confidence interval} \item{model.fit}{if ESTmethod = "SP", object of type
  'coxph' from fitting the model \code{coxph(Surv(time,
  event)~Y)}} \item{cutoff, CImethod, SEmethod,
  predict.time, alpha}{ function inputs}
}
\description{
  This function estimates the AUC, TPR(c), FPR(c), PPV(c),
  and NPV(c) for for a specific timepoint and marker cutoff
  value c using semiparametric or nonparametric estimates. 
  Standard errors, and confidence intervals are
  also computed. Either analytic or bootstrap standard
  errors can be computed.
}
\examples{
data(SimData)

#non-parametric estimates
tmp <- survAM.estimate(time =survTime, 
                       event = status, 
                       marker = Y,
                       data = SimData,
                       estimation.method = "IPW",
                       predict.time = 2, 
                       marker.cutpoint = 0, 
                       bootstraps = 50)
tmp
tmp$estimates

#semi-parametric estimates
tmp <- survAM.estimate(time =survTime, 
                       event = status, 
                       marker = Y,
                       data = SimData,
                       estimation.method = "Cox",
                       predict.time = 2, 
                       marker.cutpoint = 0, 
                       bootstraps = 50)
                       
#semi-parametric estimates with asymptotic standard errors
tmp <- survAM.estimate(time =survTime, 
                       event = status, 
                       marker = Y,
                       data = SimData,
                       estimation.method = "Cox", 
                       se.method = "asymptotic", 
                       predict.time = 2, 
                       marker.cutpoint = 0, 
                       bootstraps = 50)

}
\references{
   Liu D, Cai T, Zheng Y. Evaluating the predictive 
   value of biomarkers with stratified case-cohort 
   design. Biometrics 2012, 4: 1219-1227.


  Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance
  of markers for future events. Lifetime Data Analysis.
  2008, 14: 86-113.

  Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent
  predictive values of prognostic biomarkers with failure
  time outcome. JASA 2008, 103: 362-368.
}

