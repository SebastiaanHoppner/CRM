scaleResidualsByMAD <- function (residuals) {
  if (length(residuals) / 2 > sum(residuals == 0)) {
    residuals <- abs(residuals) / (1.4826 * median(abs(residuals)))
  } else {
    residuals <- abs(residuals) / (1.4826 * median(abs(residuals[residuals != 0])))
  } # all scaled residuals are non-negative
  return(residuals)
}
