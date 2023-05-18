AndrewsWeightFunction <- function (x, q) {
  i1 <- which(abs(x) <= pi*q)
  i2 <- which(abs(x) > pi*q)
  
  andrews_weights <- x
  andrews_weights[i1] <- (q / x[i1]) * sin(x[i1] / q)
  andrews_weights[i2] <- 0
  return(as.numeric(andrews_weights))
}