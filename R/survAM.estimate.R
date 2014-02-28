survAM.estimate <- function(time, event, marker,
                             data, 
                             predict.time,  
                             marker.cutpoint = 'median', 
                             estimation.method = "NP", 
                             ci.method = "logit.transformed", 
                             bootstraps = 1000, 
                             alpha=0.05){

  # checks
  stopifnot(is.data.frame(data))
  
  time <- eval(substitute(time), data)
  event <- 1*eval(substitute(event), data)
  marker <- eval(substitute(marker), data)

  stopifnot(is.element(estimation.method, c("NP", "SP")))
  stopifnot(is.numeric(predict.time))
  if(marker.cutpoint=='median') marker.cutpoint =  median(eval(substitute(marker), data))
  stopifnot(is.numeric(marker.cutpoint))
  
  #set some defaults
  measures = c('all')
  cutoff <- marker.cutpoint
  cutpoint <- marker.cutpoint
  cutoff.type = "none";
  SEmethod ="bootstrap"
  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")
  
  
  ESTmethod = estimation.method
  CImethod = ci.method
  N = nrow(data)
  
  
  ## get estimates via getEstimates, also calculate the bootstrap se if necessary
  
  #we handle semiparametric and nonparametric estimates differently
  if(is.element(ESTmethod, c("S", "SP", "Semi-Parametric", "semiparametric"))){
  
    mydata <- prepareDataSP(time, event, marker)
    
    if(SEmethod == "normal"){
      
      estRawOutput <- getEstimatesSP( data = mydata, 
                                      cutpoint = cutpoint,  
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = TRUE,  
                                      cutoffN = N)  
      
    }else if(substr(SEmethod, 1,4)=="boot"){
      bootstraps = round(bootstraps)
      if(bootstraps <= 1) stop("bootstraps must be larger than 1")
      #estimates
      estRawOutput<-  getEstimatesSP( data = mydata, 
                                                     cutpoint = cutpoint,  
                                                     measures = measures,
                                                     predict.time = predict.time,
                                                     CalVar = FALSE,  
                                                     cutoffN = N)
      #bootstrap ci's
      bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
      for( b in 1:bootstraps){                  
        bootests[b,] <- unlist(getEstimatesSP( data = mydata[sample.int(N, replace = TRUE),], 
                                                        cutpoint = cutpoint,  
                                                        measures = measures,
                                                        predict.time = predict.time,
                                                        CalVar = FALSE,  
                                                        cutoffN = N)$est) 
      }
      estRawOutput$se <- data.frame(t(apply(bootests, 2, sd, na.rm = TRUE)))
      names(estRawOutput$se) = names(estRawOutput$estimates)
    }
    
  }else if(is.element(ESTmethod, c("N", "NP", "Non-Parametric", "nonparametric"))){
    mydata <- prepareDataNP(time, event, marker)
    
    
    ##no bootstrap
   
      if(SEmethod == "normal"){
        
        estRawOutput <- getEstimatesNP( data = mydata, 
                                        cutpoint = cutpoint,  
                                        measures = measures,
                                        predict.time = predict.time,
                                        CalVar = TRUE,  
                                        cutoffN = N)  
        
      }else if(substr(SEmethod, 1,4)=="boot"){
        bootstraps = round(bootstraps)
        if(bootstraps <= 1) stop("bootstraps must be larger than 1")
        #estimates
        estRawOutput<-  getEstimatesNP( data = mydata, 
                                                        cutpoint = cutpoint,  
                                                        measures = measures,
                                                        predict.time = predict.time,
                                                        CalVar = FALSE  
                                                       )
        #bootstrap ci's
        bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
        for( b in 1:bootstraps){   

          bootests[b,] <-  unlist(getEstimatesNP( data = mydata[sample.int(N, replace = TRUE),], 
                                                          cutpoint = cutpoint,  
                                                          measures = measures,
                                                          predict.time = predict.time,
                                                          CalVar = FALSE)$est) 
        }
 
        estRawOutput$se <- data.frame(t(apply(bootests, 2, sd, na.rm = TRUE)))
        names(estRawOutput$se) = names(estRawOutput$estimates)
      }
    
      
    
    
  }else{
    
    stop("ESTmethod not set correctly: it must be one of `SP` (or 'semiparametric') or 'NP' (or 'nonparametric')")
  }
     

     #process the raw estimate data for clean output
     myests <- processRawOutput(estRawOutput, CImethod, alpha)
 
        
  myests$cutpoint = cutpoint; 
  myests$ESTmethod = ESTmethod; 
  myests$CImethod = CImethod; 
  myests$SEmethod = SEmethod;
  myests$predict.time = predict.time; 
  myests$alpha = alpha; 
  
  ## define class
  class(myests) <-  "SurvAM"
  myests

}
