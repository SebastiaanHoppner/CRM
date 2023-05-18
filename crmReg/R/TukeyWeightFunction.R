TukeyWeightFunction <- function (x, q) {
  i1 <- which(abs(x) <= q)
  i2 <- which(abs(x) > q)
  
  tukey_weights <- x
  tukey_weights[i1] <- (1 - (x[i1] / q)^2)^2
  tukey_weights[i2] <- 0
  return(as.numeric(tukey_weights))
}