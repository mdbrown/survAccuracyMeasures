

WGT.FUN <- function(newdata, data, w.ptb=NULL, t0)
{
  ## ====================================##
  ## KM Estimator of Censoring Survival  ##
  ## ====================================##
  Ghat.FUN <- function(tt, Ti, Di,type='fl',w.ptb=NULL)
  {
    tmpind <- rank(tt); if(is.null(w.ptb)){w.ptb=rep(1,length(Ti))}
    summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type=type, weight=w.ptb), sort(tt))$surv[tmpind]
  }
  Ti = data[,1]; Di = data[,2]; tj = newdata[,1]; Wj = dj = newdata[,2]
  Wj[tj<=t0] = dj[tj<=t0]/Ghat.FUN(tj[tj<=t0],Ti,Di)
  Wj[tj >t0] = 1/Ghat.FUN(t0,Ti,Di)
  Wj
}


