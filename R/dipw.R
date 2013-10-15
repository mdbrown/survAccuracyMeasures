
#########################functions 




## ===============================================================##
## 			 Double IPW Estimate of the ROC Curve  		            ##
## ===============================================================##
EstROC.DIPW.NP.FUN <- function(data0,u0,type,c0=NULL,rtn="EST")
{
## data0: cohortdata cbind(xi,di,yi,vi,zi,si,psi,wi) with si  = stratum by di and zi; psi: marginal stratum probability 
## set vi, zi, si, psi, wi all to 1 for a full cohort analysis 
	  N = nrow(data0); ## cohort size
   	  data = data0[data0$vi==1,]  ## sampled data
     data = data [order(data$yi),] ## sorted by marker
      ck = data$yi; 
      wgtk = data$wi;  
	   xk = data$xi; 
      sk = data$si; 
      psk=data$psi; 
 
	   nc = length(ck); # sampled size 
	   ind.ck = (1:nc)[order(ck)]
      CWk = WGT.FUN(data[,c(1,2)],data0)  # use full cohort xi di to calculate censoring weights

    	St0.Fck = sum.I(ck,">=",ck,wgtk*CWk*(xk >= t0))/sum(CWk*wgtk)
    	Ft0.Fck = sum.I(ck,">=",ck,wgtk*CWk*(xk <  t0))/sum(CWk*wgtk)
    	Fck = sum.I(ck,">=",ck,wgtk*CWk)/sum(CWk*wgtk)
    	St0 = max(St0.Fck)            ## St0     = P(T> t0)
   	   Ft0 = 1-St0                   ## Ft0     = P(T<=t0)
    	FPR.ck= (St0-St0.Fck)/St0     ## P(Y> ck|T> t0)
    	TPR.ck= (Ft0-Ft0.Fck)/Ft0     ## P(Y> ck|T<=t0)
    	NPV.ck= St0.Fck/Fck           ## P(T> t0|Y<=ck)
    	PPV.ck= (Ft0-Ft0.Fck)/(1-Fck) ## P(T<=t0|Y> ck)
   	    AUC = sum(TPR.ck*(FPR.ck-c(FPR.ck[-1],0)))
	
       ## acc.ck: accuracy estimates at c
    	nm.acc = c("FPR","TPR","NPV","PPV"); 
    	acc.ck = data.frame("cutoff"=ck,"FPR"=FPR.ck,"TPR"=TPR.ck,"NPV"=NPV.ck, "PPV"=PPV.ck)    
    	
    	## acc.uk: accuracy estimate at u
    	if (!is.null(u0)) { 
    		acc.uk = acc.ck; ind0 = match(type,names(acc.ck)); 
    		uk = acc.uk[,ind0]; acc.uk = acc.uk[order(uk),]; uk = sort(uk); 
    		if(ind0==1){tmpind.u0 = sum.I(u0,">=",uk)}else{tmpind.u0 = sum.I(u0,">",uk)}
    	}
	    acc.c0=F.c0=NULL
	    if(!is.null(c0)){
    		tmpind.c0 = sum.I(c0, ">=", ck); acc.c0 = acc.ck[tmpind.c0,]; F.c0 = Fck[tmpind.c0]
    	}  #else{ acc.c0 = acc.ck }
    	est= list("AUC" = AUC, "ACC.u0"=acc.uk[tmpind.u0,-c(1,ind0)],"ACC.c0"=acc.c0)
    	#est= list("AUC" = AUC, "ACC.u0"=acc.uk[tmpind.u0,-c(1,ind0)],"ACC.c0"=acc.c0, "ACC.all" = acc.ck) ##output all cutoff
    	if(rtn=="EST"){
    		return(est)
    	}else{
    ###### Variance calculation below ########
    Phi=Phi.C.new.FUN(xk=data$xi,dk=data$di, Ti=data0$xi, Di=data0$di)
    
    ## doing u0 and c0 together
    if (!is.null(u0)) {
		c.u0 = acc.uk[tmpind.u0,1]; acc.u0.temp = acc.uk[tmpind.u0,]; 
		tmpind.u0c= sum.I(c.u0,">=",ck); F.c0.b = c(Fck[tmpind.u0c],F.c0)} else {c.u0 = NULL; acc.u0.temp=NULL; F.c0.b =F.c0} 
	 
		CC = c(c.u0,c0); nu0 = length(u0); 
		acc.c0.b = rbind(acc.u0.temp,acc.c0); 
		
	U.ACC.c0.tmp = as.list(1:4); 
	U.ACC.c0 =Wexp.c0= as.list(1:4); 
	names(U.ACC.c0.tmp) = names(U.ACC.c0)=nm.acc; n.acc.c=length(U.ACC.c0)

	
	I.ck.c0 = 1*(ck>=VTM(CC,nc)); 
    U.ACC.c0.tmp$FPR =  (xk >  t0)*(I.ck.c0-  VTM(acc.c0.b$FPR,nc))/St0      ## exp for FPRhat(c)-FPR(c)
    U.ACC.c0.tmp$TPR =  (xk <= t0)*(I.ck.c0-  VTM(acc.c0.b$TPR,nc))/(1-St0)  ## exp for TPRhat(c)-TPR(c)
	U.ACC.c0.tmp$NPV =  (1-I.ck.c0)*(1*(xk> t0)-VTM(acc.c0.b$NPV,nc))/VTM(F.c0.b,nc)
	U.ACC.c0.tmp$PPV =    I.ck.c0*(1*(xk<=t0)-VTM(acc.c0.b$PPV,nc))/(1-VTM(F.c0.b,nc))

    U.AUC = (xk<=t0)/(1-St0)*(1-FPR.ck-AUC)+(xk>t0)/St0*(TPR.ck-AUC)

	Wexp.np.AUC = CWk*U.AUC+Phi%*%(wgtk*CWk*U.AUC)/sum(wgtk)
    se.auc = sqrt(Est.Var.CCH.trueweights(N,Wexp.np.AUC,data,data$si, subcohort=TRUE))
    se.u0 = NULL
      if (!is.null(u0)) {
      	  se.u0=NULL
      	  ind1=match(type, nm.acc)
         
      	  npu = length(nm.acc)-1
      	  U.ACC.u0 = Wexp.u0 = as.list(1:npu); names(U.ACC.u0) = c(nm.acc[-ind1]); 

      	  
		  for(kk in 1:npu){
			tmpnm = names(U.ACC.u0)[kk]; 
		  	dACC.hat = dACC.FUN(u0,uu=uk,A.u = acc.uk[,match(tmpnm,names(acc.uk))], bw=NULL)
		  	U.ACC.u0[[kk]] = U.ACC.c0.tmp[[match(tmpnm,names(U.ACC.c0.tmp))]][,1:nu0] - 
					 VTM(dACC.hat,nc)*U.ACC.c0.tmp[[match(type,names(U.ACC.c0.tmp))]][,1:nu0]
              U.ACC.u0[[kk]] = CWk*U.ACC.u0[[kk]]+Phi%*%(wgtk*CWk*U.ACC.u0[[kk]])/sum(wgtk)
              se.u0 = c(se.u0,sqrt(Est.Var.CCH.trueweights(N,data.frame(U.ACC.u0[[kk]]),data,data$si, subcohort=subcohort))) 
                      
		}
		se.u0 = data.frame(matrix(se.u0,nrow=length(u0)))
   	           names(se.u0) = c(nm.acc[-ind1]) 
      }
      se.c0 = NULL
	if (!is.null(c0)) {
		npc = length(U.ACC.c0)	
		for(kk in 1:npc){ 
			if (!is.null(u0)) {U.ACC.c0[[kk]] = U.ACC.c0.tmp[[kk]][,-(1:nu0)]} else {
			U.ACC.c0[[kk]] = U.ACC.c0.tmp[[kk]]}
       		Wexp.c0[[kk]] = CWk*U.ACC.c0[[kk]]+Phi%*%(wgtk*CWk*U.ACC.c0[[kk]])/sum(wgtk)
       		se.c0 = c(se.c0,sqrt(Est.Var.CCH.trueweights(N,data.frame(Wexp.c0[[kk]]),data,data$si, subcohort=subcohort)))
    	}
    	 se.c0 = data.frame(matrix(se.c0,nrow=length(c0)))
   	    names(se.c0) = nm.acc    
    }
    #jjunk = list(Wexp.np.AUC = Wexp.np.AUC, Wexp.np.c0=Wexp.c0, Wexp.np.u0 = Wexp.u0) 
   	 	   #Wexp = data.frame(jjunk$Wexp.AUC,jjunk$Wexp.np.c0,jjunk$Wexp.np.u0)
   	 	   #Wexp = data.frame(jjunk$Wexp.np.c0[[1]])
   	   
       list(estimates = est,"se.auc" =se.auc,"se.u"=se.u0, "se.c"=se.c0) 
       }	
}



Phi.C.new.FUN<-function(xk,dk,Ti, Di)
{
   #xk=xi; #dk=data[,7]; 
   tt=pmin(xk,t0); 
   TT=sort(unique(pmin(Ti[Di==0], t0)));
   nk=length(xk); NN=length(Ti)
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

##----------------------------------------------------------------------------------------------
## Simulate time-to-event data from											
## T: \lambda(t|y)=\lambda_0(t)exp(\beta*y)									
##    \lambda_0(t)=0.1													
##    \beta	=log(3)													
## Y: marker, N(0,1)												
## Z: Binary variable, (aux.z>0) (aux.z<=0), where aux.z~N(0,1) and Cov(Y,aux.z)=0.5		
## C: indep Cen: min( Gamma(2.5,2),2)												
##    Dep Cen:   min (exp(y/10)*Gamma(2.5, 2). 2)  
##----------------------------------------------------------------------------------------------

SIM.CCH.Z.FUN <- function(nn=5000, DepCen=F, ncch0.z0=NULL,ncch0.z1=NULL,ncch1.z0=NULL,ncch1.z1=NULL)
{
    	junk = cbind(rnorm(nn,0,1),rnorm(nn,0,1))%*%chol(matrix(c(1,0.5,0.5,1),2,2))
    	yi = junk[,1]
	yi = (yi>=0)*pmin(abs(yi), 5)-(yi<0)*pmin(abs(yi), 5)
#    	zi = ifelse(junk[,2]>0,1,0) #oldstra
    	zi = ifelse(junk[,1]>0,1,0) # newstra
#    	beta = log(3)
    	ti = -beta*yi+log(-log(runif(nn)));ti<-exp(ti)*10
    	if(!DepCen){
        ci = pmin(rgamma(nn,2.5,2),2)
    	}else{
        ci = pmin(exp(yi/10)*rgamma(nn,2,2),2)
    	}
 	xi = pmin(ti,ci); di = 1*(ti<=ci);  vi =rep(0,nn); si = rep(0, nn); psi = rep(0, nn);

	ncch0=ncch0.z0+ncch0.z1
	ncch1=ncch1.z0+ncch1.z1

	case.z1.ind = (1:nn)[di==1 & zi==1]
      case.z0.ind = (1:nn)[di==1 & zi==0]
      control.z1.ind = (1:nn)[di==0 & zi==1]
    	control.z0.ind = (1:nn)[di==0 & zi==0]
     	ncase.z1 = length(case.z1.ind)
      ncase.z0 = length(case.z0.ind)
      ncontrol.z1 = length(control.z1.ind)
    	ncontrol.z0 = length(control.z0.ind)
      temp.ncch10 = min(ncase.z0,ncch1.z0)
      temp.ncch11 = min(ncase.z1,ncch1.z1)  
      temp.ncch00 = min(ncontrol.z0,ncch0.z0+(ncch1.z0-temp.ncch10))
      temp.ncch01 = min(ncontrol.z1,ncch0.z1+(ncch1.z1-temp.ncch11))

    	junk.ind00 = sample(control.z0.ind, temp.ncch00,replace=F)
    	vi[junk.ind00] = 1
      junk.ind01 = sample(control.z1.ind, temp.ncch01,replace=F)
    	vi[junk.ind01] = 1 
      junk.ind10 = sample(case.z0.ind, temp.ncch10,replace=F)
    	vi[junk.ind10] = 1
    	junk.ind11 = sample(case.z1.ind, temp.ncch11,replace=F)
    	vi[junk.ind11] = 1 
	
	si[control.z0.ind] =1; 
      si[control.z1.ind] =2;  
      si[case.z0.ind] = 3; 
      si[case.z1.ind] = 4; 

	psi[control.z0.ind] =ncontrol.z0/nn; 
      psi[control.z1.ind] =ncontrol.z1/nn; 
      psi[case.z0.ind] = ncase.z0/nn; 
      psi[case.z1.ind] = ncase.z1/nn; 

	cohortdata=data.frame(xi=xi,di=di,yi=yi,vi=vi,zi=zi,si=si,psi=psi)
    	cohortdata
}
##----------------------------------------------------------------------------------------------##
## Obtain samplying probablity	`											
##----------------------------------------------------------------------------------------------##

P0HAT.cch.z.FUN.type2<-function(data) {
	xi = data[,1]; di = data[,2]; vi = data[,4]; phati = vi; zi = data[,5]
	tmp.ncch1.z1=sum((vi==1)*(zi==1)*(di==1))
	tmp.ncch1.z0=sum((vi==1)*(zi==0)*(di==1))
	tmp.ncch0.z1=sum((vi==1)*(zi==1)*(di==0))
	tmp.ncch0.z0=sum((vi==1)*(zi==0)*(di==0))

	p0hat =  di*(zi*tmp.ncch1.z1/sum(zi*di) + (1-zi)*tmp.ncch1.z0/sum((1-zi)*di))+
              (1-di)*(zi*tmp.ncch0.z1/sum(zi*(1-di)) + (1-zi)*tmp.ncch0.z0/sum((1-zi)*(1-di)))
	phati = vi*p0hat
}

##----------------------------------------------------------------------------------------------##
## Calculate Robust variance matrix for N iid vector U_ij from subject i in stratum j			
##----------------------------------------------------------------------------------------------##


WGT.FUN <- function(newdata, data,w.ptb=NULL)
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

Kern.FUN <- function(zz,zi,bw,kern0="gauss") ## returns an (n x nz) matrix ##
{ 
  out = (VTM(zz,length(zi))- zi)/bw
  switch(kern0,
         "epan"= 0.75*(1-out^2)*(abs(out)<=1)/bw,
         "gauss"= dnorm(out)/bw)
}

##------------------------------------------------------------
## calculates the derivative function of the accuracy measure 
##------------------------------------------------------------

dACC.FUN <- function(u0, uu=SE.yy, A.u = Sp.yy, bw=NULL)
{
  data = cbind(uu,A.u); data=na.omit(data); data = data[rowSums(abs(data))<Inf,]
  uu=data[,1]; A.u=data[,2]; n.u = length(uu)
  A.u = A.u[order(uu)]; uu = sort(uu)
  if(is.null(bw)){bw=1.06*min(sd(uu),IQR(uu)/1.34)*n.u^(-bw.power)}
  Ki.u0 = Kern.FUN(u0, uu[-1], bw) ## (n.uu-1) x n.u0
  c(t(A.u[-1]-A.u[-n.u])%*%Ki.u0)
}

