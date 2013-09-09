EST.Sy<-function(data, predict.time) {
  
  Y.old  <- as.matrix(data[,!is.element(names(data), c("di", "xi", "wi","vi"))])
  
  data <- data[order(-data$xi),]
  Y  <- as.matrix(data[,!is.element(names(data), c("di", "xi", "wi","vi"))])
  
  fit<-coxph( Surv(data$xi,data$di) ~ Y, weights = data$wi)
  
  beta<-fit$coef
  linearY <- Y%*%beta 
  
  r.riskset <-data$wi/cumsum(exp(linearY)*data$wi)
  
  Lambda0t  <- sum(r.riskset[(data$xi <= predict.time)&(data$di == 1)])
  
  linearY <- Y.old%*%beta
  Sy <- exp(-Lambda0t*exp(linearY))
  
  list(beta,Sy, linearY)
}
