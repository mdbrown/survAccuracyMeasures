survAM.estimate <- function(time, event, marker, predict.time, measures = c('all'), 
                   cutpoint = median(marker), CImethod = "logit.transformed", 
                   SEmethod ="normal", bootstraps = 1000, alpha=0.05){
  cutoff <- cutpoint
  cutoff.type = "none"; cutoffN = 100 #functionality to be added later
  weights = NULL
  #put checks here
  if(length(cutoff)==0) cutoff = NA;  
  
  #CImethod is either "standard" or "logit.transformed" 
  if(!is.element(substr(CImethod, 1,4), c("stan", "logi"))) stop("CImethod must be either 'standard' or 'logit.transformed'")

  #SEmethod is either "normal" or "boostrap"
  if(!is.element(substr(SEmethod, 1,4), c("norm", "boot"))) stop("SEmethod must be either 'normal' or 'bootstrap'")
  
  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")
  
  #make sure we have a cutoff if the measures call for it
  if(any(c("TPR", "FPR", "PPV", "NPV") %in% measures) & is.na(cutoff)) stop("'cutoff' must be set in order to calculate 'FPR', 'TPR', 'PPV, 'NPV'")
  N <- dim(time)[1]; if(is.null(N)) N = length(time)
  
  if(is.null(weights)){ 
    weights = rep(1, N)
    subcohort = FALSE
  }else{
    if(any(weights <=0)) stop("weights must be > 0.")
    if(length(weights)!=N) stop("length of weights must be the same as the length of time, event and marker")
    if(is.element(substr(SEmethod, 1,4), c("boot"))) stop("bootstrap SE's cannot be calculated when sample weights are provided, please set SEmethod='normal'")
    subcohort = TRUE
  }
  
  #end of checks
  
  
  ## get estimates via getEstimates, also calculate the bootstrap se if necessary
  
  #build data frame for getEstimates
  mydata <- as.data.frame(cbind(time, event, marker))
  names(mydata) = c("xi", "di", "Y")
  mydata$wi = weights
  
  if(SEmethod == "normal"){

    myests <- getEstimates( data = mydata, cutpoint = cutoff,  measures = measures, predict.time = predict.time, CalVar = TRUE, cutoff.type = cutoff.type, cutoffN, subcohort)
   
  }else if(substr(SEmethod, 1,4)=="boot"){
    bootstraps = round(bootstraps)
    if(bootstraps <= 1) stop("bootstraps must be larger than 1")
    myests <- getEstimates( data = mydata, cutpoint = cutoff,  measures = measures, predict.time = predict.time, CalVar = FALSE)
    
    #bootstrap ci's
    bootests <- matrix(ncol = length(myests$est), nrow = bootstraps)
    for( b in 1:bootstraps){
      bootests[b,] <- getEstimates( data = mydata[sample(1:N, replace = TRUE),], cutpoint = cutoff,  measures = measures, predict.time = predict.time, CalVar = FALSE)$est    
    }
    myests$se <- apply(bootests, 2, sd)
  }
  
   
  #calculate confidence intervals
  if(substr(CImethod, 1, 4)=="stan"){
    
    myests$CIbounds = data.frame(rbind(upper = myests$estimates - qnorm(alpha/2)*myests$se, 
                 lower = myests$estimates - qnorm(1-alpha/2)*myests$se))
    names(myests$CIbounds) = c("coef", measures)
    
  }else{
    #logit transform everything but the coef
  
    myests$CIbounds = data.frame(rbind(upper = expit(logit(myests$estimates[-1]) - qnorm(alpha/2)*myests$se[-1]/(myests$estimates[-1]*(1-myests$estimates[-1]))), 
                                       lower = expit(logit(myests$estimates[-1]) - qnorm(1-alpha/2)*myests$se[-1]/(myests$estimates[-1]*(1-myests$estimates[-1])))))
    myests$CIbounds = cbind(data.frame(rbind(upper = myests$estimates[1] - qnorm(alpha/2)*myests$se[1], 
                                       lower = myests$estimates[1] - qnorm(1-alpha/2)*myests$se[1])), myests$CIbounds)
    names(myests$CIbounds) = c("coef", measures)
    
  }
  names(myests$estimates) = c("coef", measures)
  names(myests$se) = c("coef", measures)
  myests$model.fit <- myests$fit; 
  myests$cutpoint = cutoff; 
  myests$CImethod = CImethod; 
  myests$SEmethod = SEmethod;
  myests$predict.time = predict.time; 
  myests$alpha = alpha; 
  
  ## return the results in a nice fashion
  class(myests) <-  "SurvAM"
  myests
}
