##### SEMI-PARAMETRIC

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

#deprecated
Est.Wexp.cpp <-
  function(data,N,RT.out,predict.time,uu0Vec,typexVec,typeyVec, resid.sco, fit.var, cutoffs) {
    
    if(missing(data))      { stop("Est.Wexp0: data not specified") }  
    
    #if( !("status" %in% names(data)) )  { stop(sprintf(errFormat,"status")) }
    #if( !("times" %in% names(data)) )  { stop(sprintf(errFormat,"times")) }
    
    numCuts = nrow(data)
    nr = numCuts
    if(!"wi" %in% names(data)) {
      data$weights=1
    }
    
    # First, fit the survival model    
    data = data[order(data$linearY),]   ## it is sorted it before this function; 
    
    Y  <- as.matrix(data[,!is.element(names(data), c("times", "zi", "status", "wi", "vi","Sy","linearY"))])
    
    np = dim(Y)[2]
    # fit  = coxph(Surv(data$times,data$status)~Y, 
    #              method="breslow", weight=data$weights)   
    
    # Doing riskmat, haz0 and time by hand since coxph.detail appears 
    #  to be a newer R feature & some users may have not updated their R.
    #    Note: this hazard is frequently normalized,
    #    by multiplying by exp(mean(data$Y)*fit$coef), but that is 
    #    not necessary here, as our haz0 below doesn't want it.
    
    
    #dataD    =  subset(data[order(data$times),],status==1)  

    if(is.na(cutoffs)[1]){
    Wexp.all <- getWEXP(as.matrix(data), as.matrix(Y), N, as.matrix(RT.out), predict.time, c(resid.sco), fit.var);
    }else{
      
      cutpos = sum.I(cutoffs,">=", data$linearY)
      
      Y.sub <- as.matrix(Y[cutpos,])
      subdata <- data[cutpos,]
     
      ncut = nrow(subdata)
      
      np = dim(Y)[2]
      

    Wexp.all <- getWEXPcutoff(as.matrix(data),
                               as.matrix(subdata),
                               Y = as.matrix(Y),
                               Y.sub,
                               N, as.matrix(RT.out), 
                               predict.time, c(resid.sco), fit.var, 
                               cutoffs);

      
    }
    ## now get iid expansion for other accuracy summaries
    ## global summaries 
    ## AUC = sum(RT.out$TPR*(RT.out$FPR-c(RT.out$FPR[-1],0)))
    mmm = length(RT.out$TPR)
    #ITPR = sum(RT.out$TPR*(RT.out$RiskT-c(0,RT.out$RiskT[-mmm])))
    #IFPR = sum(RT.out$FPR*(RT.out$RiskT-c(0,RT.out$RiskT[-mmm])))
    #IDI = ITPR - IFPR
    
    #Wexp.ITPR = Wexp.all[[4]]%*%(RT.out$RiskT-c(0,RT.out$RiskT[-mmm]))+
    #             (Wexp.all[[1]]-cbind(0,Wexp.all[[1]][,-mmm]))%*%RT.out$TPR
    #Wexp.IFPR = Wexp.all[[3]]%*%(RT.out$RiskT-c(0,RT.out$RiskT[-mmm]))+
    #             (Wexp.all[[1]]-cbind(0,Wexp.all[[1]][,-mmm]))%*%RT.out$FPR
    #Wexp.IDI= Wexp.ITPR - Wexp.IFPR 	
    #Wexp.AUC = Wexp.all[[4]]%*%(RT.out$FPR-c(RT.out$FPR[-1],0))+(Wexp.all[[3]]-cbind(Wexp.all[[3]][,-1],0))%*%RT.out$TPR
    Wexp.AUC = cbind(0,Wexp.all[[4]])%*%(c(1,RT.out$FPR)-c(RT.out$FPR,0))+
      (cbind(0,Wexp.all[[3]])-cbind(Wexp.all[[3]],0))%*%c(1,RT.out$TPR)
    
    
    
    
    if(!is.null(uu0Vec)){
      nvp = length(uu0Vec)
      Wexp.vp  = matrix(0,nr,nvp)
      
      for (pp in 1:nvp) {	
        uu0 = uu0Vec[pp]   
        
        uuk = sort(RT.out[,1]); 
        tmpind = sum.I(uu0,">=",uuk)
        ind0.y = match(typeyVec[pp],c("RiskT","v","FPR","TPR","rho","NPV","PPV"))
        
        Wexp.vp[,pp] = Wexp.all[[ind0.y]][,tmpind]
      }
      
    }else{
      Wexp.vp = NULL
    }
    
    
    list(Wexp.beta = Wexp.all[[8]], Wexp.AUC = Wexp.AUC,Wexp.vp=Wexp.vp)   
  }


Est.Wexp<-function(data,N,RT.out,predict.time,uu0Vec,typexVec,typeyVec, resid.sco, fit.var) {
  
  if(missing(data))      { stop("Est.Wexp0: data not specified") }  
  
  # if( !("status" %in% names(data)) )  { stop(sprintf(errFormat,"status")) }
  # if( !("times" %in% names(data)) )  { stop(sprintf(errFormat,"times")) }
  
  numCuts = nrow(data)
  nr = numCuts
  
  
  # First, fit the survival model    
  data = data[order(data$linearY),]   ## it is sorted it before this function; 
  #data = data[order(data$Sy),] 
  Y  <- as.matrix(data[,!is.element(names(data), c("times", "zi", "status", "wi", "vi","Sy","linearY"))])
  
  np = dim(Y)[2]
  
  #fit  = coxph(Surv(data$times,data$status)~Y, 
  #            method="breslow", weights=data$wi)   
  
  # Doing riskmat, haz0 and time by hand since coxph.detail appears 
  #  to be a newer R feature & some users may have not updated their R.
  #    Note: this hazard is frequently normalized,
  #    by multiplying by exp(mean(data$Y)*fit$coef), but that is 
  #    not necessary here, as our haz0 below doesn't want it.
  status <- NULL
  rrk      =  exp(data$linearY)
  dataD    =  subset(data[order(data$times),],status==1)
  riskmat  =  t(sapply(data$times,function(x) x >= dataD$times))
  
  s0   = t(riskmat) %*% (rrk*data$wi) ## length of nt
  s1   = t(riskmat) %*% t(VTM(rrk*data$wi,np)*t(Y))  ## nt *np 
  haz0      = dataD$wi / colSums(riskmat*rrk*data$wi)
  cumhaz0   = cumsum(haz0)
  cumhaz.t0 = sum.I(predict.time, ">=", dataD$times, haz0)
  ## CondSyk   = exp(-cumhaz.t0*rrk) ## check it is the same as Sy 
  tmpind    = (data$times<=predict.time)&(data$status==1)
  tmpind.t  = sum.I(data$times[tmpind], ">=", dataD$times)
  
  #resid.sco         = resid(fit, type="score")
  Wexp.beta         = resid.sco %*% fit.var * N
  Wexp.Lam1         = rep(0, numCuts)
  Wexp.Lam1[tmpind] = N/s0[tmpind.t]
  Wexp.Lam1 = Wexp.Lam1 - sum.I(pmin(predict.time,data$times), ">=", dataD$times, haz0/s0)*rrk*N
  Wexp.Lam2 = Wexp.beta %*% sum.I(predict.time, ">=", dataD$times, haz0*s1/t(VTM(s0,np)))
  Wexp.Lam  = Wexp.Lam1 - Wexp.Lam2
  
  # end of most basic expansions... 
  #    next calcs are derived expansions of performance measures
  # Fyk  = sum.I(data$Sy, ">=", data$Sy, data$wi)/sum(data$wi)  ##Fyk = P(Sy <c)
  Fyk  = sum.I(data$linearY, ">=", data$linearY, data$wi)/sum(data$wi)  # Fyk is the distribution of linear predictor under the cox model, same if use -Sy but not Sy
  #Fyk   = rank(data$Y,ties="max")/numCuts
  dFyk = Fyk - c(0,Fyk[-numCuts])
  
  St0.Fyk   = cumsum(data$Sy*dFyk)  ## St0.Fyk = P(T> t0,Sy<=c)
  St0       = max(St0.Fyk)          ## St0 = P(T>t0)
  St0.Syk   = St0-St0.Fyk           ## St0.Syk = P(T>t0,Sy>c) 
  
  Wexp.Cond.Stc = -VTM(data$Sy*rrk, numCuts) *
    (t(VTM(Wexp.Lam,numCuts))+cumhaz.t0*Wexp.beta%*%t(Y)) ## iid expansion of Sy at each c=sy;
  ## changed below to >= from < 
  Wexp.Stc  = t(sum.I(data$linearY, "<", data$linearY, t(Wexp.Cond.Stc)*dFyk)) + 
    data$Sy*(data$linearY > VTM(data$linearY,numCuts)) - 
    VTM(St0.Syk, numCuts)     ## iid expansion of St0.Syk  
  
  Wexp.St   = colSums(t(Wexp.Cond.Stc)*dFyk) + data$Sy - St0  ## iid expansion of St0, 
  
  Wexp.Fc   = 1*(data$linearY <= VTM(data$linearY,numCuts)) - VTM(Fyk,numCuts)  ## iid expansion of Fyc = P(Sty<c) for c known;  
  
  ## Assemble for classic performance measures: given linear predictor; 
  Wexp.all = as.list(1:7); 
  
  names(Wexp.all)=c("RiskT","v","FPR","TPR","rho","NPV","PPV")
  
  Wexp.all[[1]] = -Wexp.Cond.Stc; ## iid expansion of conditional risk: Fty = P(T<t|Y) for each person n*n
  Wexp.all[[2]] = Wexp.Fc; ## iid expansion of v=Fyc = P(linearY<c)
  Wexp.all[[3]] = ( - Wexp.St * VTM(RT.out[,4], numCuts) + Wexp.Stc)/St0  ## iid expansion of P(Sy>c|T>t)
  Wexp.all[[4]] = (Wexp.St* VTM(RT.out[,5], numCuts) -Wexp.Fc - Wexp.Stc)/(1-St0)  ## iid expansion of P(Sy>c|T<t)
  Wexp.all[[5]] =  -Wexp.St;  ## iid expansion of P(T>t)
  Wexp.all[[6]] = (Wexp.St-Wexp.Stc-VTM(RT.out[,7], nr)*Wexp.Fc)/VTM(Fyk,nr)                
  Wexp.all[[7]]= (VTM(RT.out[,6]-1, nr)*Wexp.Fc-Wexp.Stc)/VTM(1-Fyk,nr)
  ## now get iid expansion for other accuracy summaries
  ## global summaries 
  ## AUC = sum(RT.out$TPR*(RT.out$FPR-c(RT.out$FPR[-1],0)))
  mmm = length(RT.out$TPR)
  
  Wexp.AUC = cbind(0,Wexp.all[[4]])%*%(c(1,RT.out$FPR)-c(RT.out$FPR,0))+
    (cbind(0,Wexp.all[[3]])-cbind(Wexp.all[[3]],0))%*%c(1,RT.out$TPR)
  
  
  if(!is.null(uu0Vec)){
    nvp = length(uu0Vec)
    Wexp.vp  = matrix(0,nr,nvp)
    
    for (pp in 1:nvp) {	
      uu0 = uu0Vec[pp]   
      
      uuk = sort(RT.out[,1]); 
      tmpind = sum.I(uu0,">=",uuk)
      ind0.y = match(typeyVec[pp],c("RiskT","v","FPR","TPR","rho","NPV","PPV"))
      
      Wexp.vp[,pp] = Wexp.all[[ind0.y]][,tmpind]
    }
    
  }else{
    Wexp.vp = NULL
  }
  
  
  
  list(Wexp.beta = Wexp.beta, Wexp.AUC = Wexp.AUC, Wexp.vp=Wexp.vp)   
}

#this goes in dipw.R




##non parametric


Phi.C.new.FUN<-function(xk,dk,Ti, Di, t0){
      #xk=xi; #dk=data[,7]; 
        
        tt=pmin(xk,t0); 
        TT=sort(unique(pmin(Ti[Di==0], t0)));
        nk=length(xk); N=length(Ti)
        junk=summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl'), TT)
        pi=junk$n.risk/N
        dLambda=junk$n.event/junk$n.risk
        #c(0, diff(junk$surv))
          tmp.ti=rep(xk, each=nk)
          tmp.tj=rep(xk, nk)
          tmp.t=pmin(tmp.ti, tmp.tj, t0)
          phi2=matrix(sum.I(tmp.t, ">=", TT, dLambda/pi), nk, nk)
          tmpind <- rank(tt); 
          pk=summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type='fl'), sort(tt))$n.risk[tmpind]/N
          phi1=matrix((tmp.ti<=tmp.tj)*rep(1-dk, each=nk)/rep(pk, nk), nk, nk)
          phi=phi1-phi2
          t(phi)
          #  row of the output is for subject
            #  colum of the output is for time
          }


