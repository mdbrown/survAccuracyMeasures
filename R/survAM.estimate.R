survAM.estimate <- function(time, event, marker,
                             data, 
                             predict.time,  
                             threshold, 
                             threshold.type = c("marker", "FPR", "TPR", "PPV", "NPV"),
                             estimation.method = c("IPW", "Cox"), 
                             ci.method = c("logit.transformed", "standard"), 
                             se.method = c("bootstrap", "asymptotic"),
                             bootstraps = 1000, 
                             alpha=0.05){

  # checks
  stopifnot(is.data.frame(data))
  
  time <- eval(substitute(time), data)
  event <- 1*eval(substitute(event), data)
  marker <- eval(substitute(marker), data)
  cutpoint.type <- threshold.type
  estimation.method <- match.arg(estimation.method)
  stopifnot(is.numeric(predict.time))
  threshold.type <- match.arg(threshold.type)
  ci.method <- match.arg(ci.method)
  stopifnot(is.numeric(threshold))
  se.method <- match.arg(se.method)
  if(threshold.type != "marker" ) {
    if(threshold <= 0 | threshold >= 1 ) stop("threshold must be between 0 and 1 if threshold.type is not 'marker'")
  }
  if(threshold.type != "marker" & se.method == "asymptotic"){
    stop("asymptotic se calculations only available for threshold.type = 'marker' ")
    
  }

  #cant return IPW se estimates 
 # if(estimation.method =="IPW" & se.method=="asymptotic") stop("Asymptotic variance calculations are not available for IPW estimates, please use the bootstrap to calculate standard error")
  
  #set some defaults
  measures = c('all')

  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")

  N = nrow(data)

  ## get estimates via getEstimates, also calculate the bootstrap se if necessary
  
  #we handle semiparametric and nonparametric estimates differently
  if(is.element(estimation.method, c("S", "Cox", "Semi-Parametric", "semiparametric"))){
  
    mydata <- prepareDataSP(time, event, marker)
    
    if(se.method == "asymptotic"){
      
      estRawOutput <- getEstimatesSP( data = mydata, 
                                      cutpoint = threshold,  
                                      cutpoint.type = cutpoint.type, 
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = TRUE)  
      
    }else if(substr(se.method, 1,4)=="boot"){
      bootstraps = round(bootstraps)
      if(bootstraps <= 1) stop("bootstraps must be larger than 1")
      #estimates
      estRawOutput<-  getEstimatesSP( data = mydata, 
                                      cutpoint = threshold,
                                      cutpoint.type = cutpoint.type, 
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = FALSE)
      #bootstrap ci's
      bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
      for( b in 1:bootstraps){                  
        bootests[b,] <- unlist(getEstimatesSP( data = mydata[sample.int(N, replace = TRUE),], 
                                               cutpoint = threshold,
                                               cutpoint.type = cutpoint.type,   
                                               measures = measures,
                                               predict.time = predict.time,
                                               CalVar = FALSE)$est) 
      }
      estRawOutput$se <- data.frame(t(apply(bootests, 2, sd, na.rm = TRUE)))
      names(estRawOutput$se) = names(estRawOutput$estimates)
    }
    
  }else if(is.element(estimation.method, c("NP", "IPW"))){
    mydata <- prepareDataNP(time, event, marker)
    
    
    ##no bootstrap
   
      if(se.method == "asymptotic"){
        
        estRawOutput <- getEstimatesNP( data = mydata, 
                                        cutpoint = threshold,  
                                        cutpoint.type = cutpoint.type, 
                                        measures = measures,
                                        predict.time = predict.time,
                                        CalVar = TRUE)  
        
      }else if(substr(se.method, 1,4)=="boot"){
        bootstraps = round(bootstraps)
        if(bootstraps <= 1) stop("bootstraps must be larger than 1")
        #estimates
        estRawOutput<-  getEstimatesNP( data = mydata, 
                                                        cutpoint = threshold,  
                                        cutpoint.type = cutpoint.type, 
                                        measures = measures,
                                        predict.time = predict.time,
                                        CalVar = FALSE )
        #bootstrap ci's
        bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
        for( b in 1:bootstraps){   

          bootests[b,] <-  unlist(getEstimatesNP( data = mydata[sample.int(N, replace = TRUE),], 
                                                  cutpoint = threshold,  
                                                  cutpoint.type = cutpoint.type, 
                                                  measures = measures,
                                                  predict.time = predict.time,
                                                  CalVar = FALSE)$est) 
        }
 
        estRawOutput$se <- data.frame(t(apply(bootests, 2, sd, na.rm = TRUE)))
        names(estRawOutput$se) = names(estRawOutput$estimates)
      }
    
      
    
    
  }else{
    
    stop("estimation.method not set correctly: it must be one of `Cox` or 'IPW'")
  }
     

     #process the raw estimate data for clean output
     myests <- processRawOutput(estRawOutput, ci.method, alpha)
     myests$roc <-  estRawOutput$roc
  myests$threshold.type = threshold.type; 
  myests$threshold = threshold; 
  myests$estimation.method = estimation.method; 
  myests$ci.method = ci.method; 
  myests$se.method = se.method;
  myests$predict.time = predict.time; 
  myests$alpha = alpha; 
  
  ## define class
  class(myests) <-  "SurvAM"
  myests

}
