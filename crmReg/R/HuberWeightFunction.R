HuberWeightFunction <- function (x, q) {
  i1 <- which(abs(x) <= q)
  i2 <- which(abs(x) > q)
  
  huber_weights <- x
  huber_weights[i1] <- 1
  huber_weights[i2] <- q / abs(x[i2])
  return(as.numeric(huber_weights))
}