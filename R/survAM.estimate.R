survAM.estimate <- function(time, event, marker, predict.time, measures = c('all'), 
                   cutpoint = median(marker), ESTmethod = "NP", CImethod = "logit.transformed", 
                   SEmethod ="normal", bootstraps = 1000, alpha=0.05){
 
  cutoff.type = "none";

    
  #put checks here
  if(length(cutpoint)==0) cutpoint = NA;  
  
  #CImethod is either "standard" or "logit.transformed" 
  if(!is.element(substr(CImethod, 1,4), c("stan", "logi"))) stop("CImethod must be either 'standard' or 'logit.transformed'")

  #SEmethod is either "normal" or "boostrap"
  if(!is.element(substr(SEmethod, 1,4), c("norm", "boot"))) stop("SEmethod must be either 'normal' or 'bootstrap'")
  
  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")
  
  #make sure we have a cutpoint if the measures call for it
  if(any(c("TPR", "FPR", "PPV", "NPV") %in% measures) & is.na(cutpoint)) stop("'cutpoint' must be set in order to calculate 'FPR', 'TPR', 'PPV, 'NPV'")
  N <- dim(time)[1]; if(is.null(N)) N = length(time)
  
  if(!all.equal(c(length(time), length(event), length(marker)), rep(N, 3))) stop("time, event and marker must be numeric vectors of equal length")
  
  #end of checks
  
  
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
                                                        CalVar = FALSE,  
                                                        cutoffN = N)
        #bootstrap ci's
        bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
        for( b in 1:bootstraps){   

          bootests[b,] <-  unlist(getEstimatesNP( data = mydata[sample.int(N, replace = TRUE),], 
                                                          cutpoint = cutpoint,  
                                                          measures = measures,
                                                          predict.time = predict.time,
                                                          CalVar = FALSE,  
                                                          cutoffN = N)$est) 
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
