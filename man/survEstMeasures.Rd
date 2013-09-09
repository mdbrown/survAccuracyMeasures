\name{survEstMeasures}
\alias{survEstMeasures}
\title{Estimate accuracy measures for a risk prediction marker from survival data}
\usage{
  survEstMeasures(time, event, marker, predict.time,
    measures = c("all"), cutpoint = median(marker),
    CImethod = "logit.transformed", SEmethod = "normal",
    bootstraps = 1000, alpha = 0.05)
}
\arguments{
  \item{time}{numeric vector, time to event variable}

  \item{event}{numeric vector, indicator for the status of
  event of interest. event = 0 for censored observations,
  and event = 1 for event of interest.}

  \item{marker}{numeric vector, marker variable of
  interest}

  \item{predict.time}{numeric value of the timepoint of
  interest for which to estimate the risk measures}

  \item{measures}{character vector that can contain the
  following 'AUC', 'TPR', 'FPR', 'PPV', 'NPV' or 'all'
  (default).}

  \item{cutpoint}{numeric value indicating the value of the
  cutpoint 'c' at which to estimate 'FPR', 'TPR', 'NPV' or
  'PPV'. default is \code{median(marker)}.}

  \item{CImethod}{character string of either
  'logit.transformed' (default) or 'standard' indicating
  whether normal approximated confidence intervals should
  be calculated using logistic transformed values or the
  standard method.}

  \item{SEmethod}{character string of either
  'normal'(default) or 'boostrap' indicating whether
  standard errors should be calculated analytically or via
  the bootstrap. The method to calculate the SE
  analytically is described in the paper referenced below.}

  \item{bootstraps}{if SEmethod = 'bootstrap', number of
  bootstrap replicates to use to estimate the SE.}

  \item{alpha}{alpha value for confidence intervals.
  (1-alpha)*100% confidence intervals are provided. default
  is alpha = 0.05.}
}
\value{
  a list with components \item{estimates}{point estimates
  for risk measures} \item{se}{standard errors for
  estimates} \item{CIbounds}{bounds for (1-alpha)*100%
  confidence interval} \item{model.fit}{object of type
  'coxph' from fitting the model \code{coxph(Surv(time,
  event)~Y)}} \item{cutoff, CImethod, SEmethod,
  predict.time, alpha}{ function inputs}
}
\description{
  This function estimates the AUC, TPR(c), FPR(c), PPV(c),
  and NPV(c) for for a specific timepoint and marker cutoff
  value c. Standard errors, and confidence intervals are
  also computed. Either analytic or bootstrap standard
  errors can be computed.
}
\examples{
data(SimData)

tmp <- survEstMeasures(time =SimData$survTime, event = SimData$status, marker = SimData$Y, predict.time = 2, cutpoint = 0)
tmp

survEstMeasures(time =SimData$survTime, event = SimData$status, marker = SimData$Y, predict.time = 2, measures = c("AUC", "TPR"), SEmethod = 'bootstrap', bootstraps = 50, cutpoint = 0)
}
\references{
  Pepe MS, Zheng Y, Jin Y. Evaluating the ROC performance
  of markers for future events. Lifetime Data Analysis.
  2008, 14: 86-113.

  Zheng Y, Cai T, Pepe MS, Levy, W. Time-dependent
  predictive values of prognostic biomarkers with failure
  time outcome. JASA 2008, 103: 362-368.
}

