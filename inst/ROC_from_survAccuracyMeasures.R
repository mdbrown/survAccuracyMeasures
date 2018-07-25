## Get ROC curve from survAM.estimate
#using for loop

library(survAccuracyMeasures)
library(tidyr)
library(ggplot2)
library(dplyr)

set.seed(112233)
#simulated data for illustration
data(SimData)
head(SimData)


ROC <- tibble(cutpoint = unique(SimData$Y)) 
ROC$FPR <-  NA
ROC$TPR <- NA

for( i in 1:nrow(ROC)) {

tmp <-   survAM.estimate(time =survTime,
                  event = status,
                  marker = Y,
                  data = SimData,
                  estimation.method = "IPW",
                  se.method = "bootstrap",
                  predict.time = 1,
                  marker.cutpoint = ROC$cutpoint[i], 
                  bootstraps = 2) 

ROC[i, c("FPR", "TPR") ] = tmp$estimates[c("FPR", "TPR")]

}

ROC %>%
  arrange(FPR) %>%
  ggplot(aes(FPR, TPR)) +
  geom_step() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") + 
  theme_bw()
