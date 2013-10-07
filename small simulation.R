# simulation to check the performance of a specific example regarding 
# the power simulations for accuracy measures brought up by Anna Baron on Sept 10, 2013
# 

#need this to simulate data
source("../survAccuracyMeasuresPower/subroutines.R")

# define parameters

S.0 = .95
t.0 = 5 ; a = -log(S.0)/t.0
predict.time = 5
cens.perc = .2

# want to compare the bootstrap, analytic and empirical SE for this case
nS = 1000 # number of sims
estimates <- bootse <- analse <- numeric(nS)

for(i in 1:nS){
  #sim data
  mydata <- SIM.data.singleMarker(1000, mu = 0, Sigma = 1, beta = 0.35, lam0 = a, cens.perc = cens.perc)  
  tmp.bootests <- survAM.estimate(time = mydata$xi, event = mydata$di, marker=mydata$Y, predict.time = predict.time,
                                  SEmethod='bootstrap', measures = "AUC" , bootstraps = 1000)
  
  bootse[i] <- tmp.bootests$se[2]
  tmp.analests <- survAM.estimate(time = mydata$xi, event = mydata$di, marker=mydata$Y, predict.time = predict.time,
                                  SEmethod='normal', measures = "AUC" , bootstraps = 1000)
  
  analse[i] <- tmp.analests$se[2]
  estimates[i] <- tmp.analests$estimates[2]
  
  print(i)
}


bootse <- bootse[1:72]
estimates <- estimates[1:72]
analse <- analse[1:72]


