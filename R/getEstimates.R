#' get Risk estimates and standard errors, called from SurvRM

#'@param data data.frame of xi (surv time), di (event ind), and Y (marker)
#'@param measures character vector of measures wanted, can be a subset of c("AUC", "beta", "TPR", "FPR", "PPV", "NPV")
#'@param predict.time numeric prediction time to evaluate measures at
#'@param CalVar Logical, should standard errors be calculated?
#'@return a list consisting of 'est' with point estimates in order of 'measures' and 
#'        (if CalVar == TRUE), 'sd' with standard errors for estimates. 
 

getEstimates <- function(data, 
                       cutpoint,  
                       measures,
                       predict.time,
                       CalVar, cutoff.type = "none", cutoffN = 100)
{  
  
#  browser()
  
  N = nrow(data)
  data$vi = 1; data$wi = 1
  cutoff <- cutpoint
  
  #junk = GetRTdata(data, predict.time)  ## data.RT and data are sorted by Y 
  ###
   
    ####
    ## Est.Sy

  
     data.s <- data[order(-data$xi),] #data sorted by stime

     fit<-coxph( Surv(data$xi,data$di) ~ data$Y, weights = data$wi, method = "breslow")   #original data, sorting doesn't matter. 
     betahat<-fit$coef                                               
  
      
     r.riskset <-data.s$wi/cumsum(exp(data.s$Y*betahat )*data.s$wi)
     Lambda0t  <- sum(r.riskset[(data.s$xi <= predict.time)&(data.s$di == 1)])
  


     linearY <-  data$Y*betahat            #linearY is sorted by original data
     Sy <- exp(-Lambda0t*exp(linearY))
  
  
    ####


  ooo = order( data$Y)
  data.RT <- cbind(data$Y, Sy, data$wi)[ooo,] 
  
  Fck <- sum.I(data.RT[,1], ">=", data.RT[,1], data.RT[,3])/sum(data.RT[,3])
  
  data.RT <- cbind(data.RT[,-c(3)],Fck) 
  
 ###


 
  
  if(cutoff.type != "none"){

   # cutoffs <- unique(sort(c( cutpoint, quantile(linearY, (1:cutoffN/cutoffN), type =1, na.rm = TRUE))))
    cutoffs <- unique(sort(c( quantile(linearY, (1:cutoffN/cutoffN), type =1, na.rm = TRUE))))
    
    cutpos = sum.I(cutoffs,">=", linearY[ooo])
    subdata.RT = data.RT[cutpos, ]
    
    RT.out = EstRTall(subdata.RT) 
    
    RT.out <- RT.out[[1]]
    AUC   = sum(RT.out$TPR*(RT.out$FPR-c(RT.out$FPR[-1],0)))
    
  }else{
  
  RT.out = EstRTall(data.RT) 
  AUC       <- RT.out[[2]]  
  RT.out    <- RT.out[[1]]
  
  }
  
  
  
  if (length(measures[measures!="AUC"])>0) {
    typey = measures[measures!="AUC"]; 
    typex = rep("cutoff", length(typey))
    vp = rep(cutoff, length(typey));

    tmpind <- with(RT.out, which.min(abs(cutoff-cutpoint)))
    all.measures = RT.out[tmpind,]
    
    
    RTvp.out  = all.measures[, measures[measures!="AUC"]]
    
  }else{
    RTvp.out = NULL; typex = typey = vp = NULL; 
  }
  
  
  if("AUC" %in% measures) est = c(betahat, AUC, RTvp.out )  else est = c(betahat, RTvp.out)
  
  
  if (CalVar)  {

    subdata = cbind(data[ooo,],data.RT[,c(2)], linearY[ooo])
    names(subdata)=c("times","status","y","vi","weights","Sy","linearY")

    jjunk = Est.Wexp(subdata,N,RT.out,predict.time,vp,typex,typey, resid(fit, "score"), fit$var)
    Wexp = cbind(jjunk$Wexp.beta,jjunk$Wexp.AUC,jjunk$Wexp.vp)
    
    se = sqrt(Est.Var.CCH.trueweights(N,Wexp,subdata,subdata$status)[[1]])  
    list(estimates = unlist(est),se =se, fit = fit) 
      
		
  } else {list(estimates = unlist(est), fit = fit)}
}

