impute_outlying_cells <- function (outlier, data, outlvarz, az, wr, increment = 0) {

  ## outlier .... index
  ## data  ...... must be a data frame
  ## outlvars ... variable names in which the outlier is outlying
  ## as ......... output spadimo loadings
  ## increment .. currently not used, set to zero
  ## wr ......... caseweights

  # require(FNN)

  outlvars <- outlvarz[[as.character(outlier)]]
  as <- az[[as.character(outlier)]]
  modifiedrow <- data[outlier, ]
  if (length(outlvars) < length(as)) {
    nn_index <- get.knnx(data = as.matrix(data[wr == 1, -outlvars]),
                         query = matrix(modifiedrow[-outlvars], nrow = 1),
                         k = 2)$nn.index
  } else {
    nn_index <- get.knnx(data = as.matrix(data[wr == 1, ]),
                         query = matrix(modifiedrow, nrow = 1),
                         k = 2)$nn.index
  }
  dataclean <- data[wr == 1, ][nn_index, ]
  imputeddataclean <- colMeans(as.matrix(dataclean[, sort(outlvars)]))

  modifiedrow[outlvars + increment] <- imputeddataclean
  return(modifiedrow)
}
