HampelWeightFunction <- function (x, q1, q2, q3) {
  i1 <- which(abs(x) <= q1)
  i2 <- which(abs(x) > q1 & abs(x) <= q2)
  i3 <- which(abs(x) > q2 & abs(x) <= q3)
  i4 <- which(abs(x) > q3)

  hampel_weights <- x
  hampel_weights[i1] <- 1
  hampel_weights[i2] <- q1 / abs(x[i2])
  hampel_weights[i3] <- (q1 * (q3 - abs(x[i3]))) / ((q3 - q2) * abs(x[i3]))
  hampel_weights[i4] <- 0
  return(as.numeric(hampel_weights))
}
