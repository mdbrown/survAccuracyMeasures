


GetRTdata = function(data,predict.time) {
  Sy = EST.Sy(data, predict.time)
  beta <- Sy[[1]]
  linearY <- Sy[[3]]
  ooo = order(linearY/beta)
  data.RT <- cbind(linearY/beta, Sy[[2]], data$wi)[ooo, ]
  Fck <- sum.I(data.RT[, 1], ">=", data.RT[, 1], data.RT[, 
                                                         3])/sum(data.RT[, 3])
  data.RT <- cbind(data.RT[, -c(3)], Fck)
  list(beta = beta, data.RT = data.RT, data = data[ooo, ])
 }