crm <- function (formula, data, maxiter = 100, tolerance = 0.01, outlyingness.factor = 1,
                 spadieta = seq(0.9, 0.1, -0.1), center = "median", scale = "qn",
                 regtype = "LTS", alphaLTS = 0.5, seed = NULL, verbose = TRUE) {
  # -----------------------------------------------------------------------------------------------
  # CRM - Cellwise Robust M-regression
  # -----------------------------------------------------------------------------------------------
  # Inputs:
  #   formula               an lm-style formula object specifying which relationship to estimate.
  #   data                  the data as a data frame.
  #   maxiter   (optional)  maximum number of iterations (default is 100).
  #   tolerance (optional)  obtain optimal regression coefficients to within a certain tolerance
  #                         (default is 0.01).
  #   outlyingness.factor (optional)  numeric value, larger or equal to 1 (default).
  #                                   Only cells are altered of cases for which the original
  #                                   outlyingness (before SPADIMO) is larger than
  #                                   outlyingness.factor * oulyingness AFTER SPADIMO.
  #                                   The larger this factor, the fewer cells are imputed.
  #   spadieta  (optional)  the sparsity parameter to start internal outlying cell detection with,
  #                         must be in the range [0,1] (default is seq(0.9, 0.1, -0.1)).
  #   center    (optional)  how to center the data. A string that matches the R function
  #                         to be used for centering (default is "median").
  #   scale     (optional)  how to scale the data. Choices are "no" (no scaling) or a string
  #                         matching the R function to be used for scaling (default is "qn").
  #   regtype   (optional)  type of robust regression. Choices are "LTS" (default) or "MM".
  #   alphaLTS  (optional)  parameter used by LTS regression. The percentage (roughly) of
  #                         squared residuals whose sum will be minimized (default is 0.5).
  #   seed      (optional)  initial seed for random generator, like .Random.seed
  #   verbose   (optional)  should output be shown during the process (default is TRUE).
  #
  # -----------------------------------------------------------------------------------------------
  # Output: A list object of class "crm" containing the following elements:
  #   coefficients          a named vector of fitted coefficients.
  #   fitted.values         the fitted response values.
  #   residuals             the residuals, that is response minus fitted values.
  #   weights               the (case) weights of the residuals.
  #   data.imputed          the data as imputed by CRM.
  #   casewiseoutliers      a vector that indicates the casewise outliers with TRUE or FALSE.
  #   cellwiseoutliers      a matrix that indicates the cellwise outliers as the difference between
  #                         the original data and imputed data, both scaled and centered.
  #   terms                 the terms object used.
  #   call                  the matched call.
  #   inputs                the list of supplied input arguments.
  #   numloops              the number of iterations.
  #   time                  the number of seconds passed to execute the CRM algorithm.
  #
  # -----------------------------------------------------------------------------------------------
  #   Written by S. Serneels, BASF Corp. and S. Höppner, KU Leuven, Aug-Nov 2017.
  #   Modified by I. Ortner July (2018) and S. Höppner (2019)
  # -----------------------------------------------------------------------------------------------


  ## Start timer
  t_start <- proc.time()


  ## Original call
  call <- match.call()


  ## Check if input arguments are correct
  if (missing(formula)) {
    stop("argument 'formula' is missing, with no default")
  } else if (class(formula) != "formula") {
    stop("argument 'formula' must be a formula")
  }
  if (missing(data)) {
    stop("argument 'data' is missing, with no default")
  } else if(!is.data.frame(data)) {
    stop("argument 'data' must be a data frame")
  }
  if (any(spadieta < 0) | any(spadieta >= 1)) {
    stop("internal SPADIMO sparsity parameter must be in the range [0,1]")
  }
  if (!regtype %in% c("LTS", "MM")) {
    stop("argument 'regtype' must be either 'LTS' or 'MM'")
  }
  if (regtype == "LTS") {
    if (alphaLTS < 0.5 | alphaLTS > 1) {
      stop("argment 'alphaLTS' must be between 0.5 and 1")
    }
  } else if (regtype == "MM") {
    alphaLTS <- NULL
  }
  if (maxiter < 0) {
    stop("argument 'maxiter' must be a positive integer")
  }
  if (tolerance <= 0) {
    stop("argument 'tolerance' must be non-negative")
  }
  if (outlyingness.factor < 1) {
    stop("argument 'outlyingness.factor' must be larger or equal to 1")
  }
  if (class(verbose) != "logical") {
    stop("argument 'verbose' must be TRUE or FALSE")
  }
  inputs <- list(formula = formula, maxiter = maxiter, tolerance = tolerance,
                 outlyingness.factor = outlyingness.factor, spadieta = spadieta,
                 center = center, scale = scale, regtype = regtype, alphaLTS = alphaLTS,
                 seed = seed, verbose = verbose)


  ## Set data in correct format (i.e. matrix with response variable as first column)
  mt <- terms(formula, data = data)
  yname <- dimnames(attr(mt, "factors"))[[1]][1]
  intercept_flag <- ifelse(attr(mt, "intercept") == 1, TRUE, FALSE)
  if (intercept_flag) { # remove intercept (i.e. column of ones)
    datam <- cbind(data[[which(colnames(data) == yname)]], model.matrix(mt, data)[, -1])
  } else {
    datam <- cbind(data[[which(colnames(data) == yname)]], model.matrix(mt, data))
  }
  colnames(datam)[1] <- yname
  nRow <- nrow(datam)
  nCol <- ncol(datam)


  ## Center and scale the data
  if (verbose) {cat(" - Center and scale the data...\n")}
  datamc <- daprpr(datam, center, scale)


  ## Apply robust regression estimator
  if (verbose) {cat(paste(" - Initial", regtype, "regression...\n"))}
  if (!is.null(seed)) set.seed(seed)
  if (regtype == "LTS") {
    robreg <- ltsReg(formula, data = as.data.frame(datamc), alpha = alphaLTS)
  } else if (regtype == "MM") {
    robreg <- lmrob(formula, data = as.data.frame(datamc))
  }


  ## Get estimated regression coefficients
  intercept <- ifelse(intercept_flag, robreg$coefficients[1], 0)
  if (intercept_flag) {
    betas <- robreg$coefficients[-1]
  } else {
    betas <- robreg$coefficients
  }


  ## Calculate residuals and derive case weights with Hample weight function
  r <- robreg$residuals
  r <- scaleResidualsByMAD(r)
  wr <- HampelWeightFunction(r, q1 = qnorm(0.95), q2 = qnorm(0.975), q3 = qnorm(0.999))
  names(wr) <- 1:nRow


  ## Find cells contributing to outlyingness
  datamc_imputed <- datamc

  if (any(wr < 1)) {
    if (verbose) {cat(" - Apply SPADIMO to find outlying cells...\n")}
    outliers <- as.list(which(wr < 1)) # casewise outliers in X-space (i.e. cases with wr[i] < 1)
    suppressWarnings(
      outlvars <- lapply(outliers, function (x) {spadimo(data = datamc[, 2:nCol],
                                                         weights = wr, obs = x,
                                                         control = list(scaleFun  = Qn,
                                                                        nlatent   = 1,
                                                                        etas      = spadieta,
                                                                        csqcritv  = 0.975,
                                                                        stopearly = TRUE,
                                                                        trace     = FALSE,
                                                                        plot      = FALSE))})
    )
    outlvarz <- lapply(outlvars, function (x) {x$outlvars})
    az <- lapply(outlvars, function (x) {x$a})
    ob <- lapply(outlvars, function (x) {x$o.before})
    oa <- lapply(outlvars, function (x) {x$o.after})
    outliersx <- as.list(as.numeric(names(which(unlist(ob) > unlist(oa) * outlyingness.factor))))
    names(outliersx) <- as.character(unlist(outliersx))


    ## Update imputed data matrix
    matz <- ldply(outliersx, impute_outlying_cells,
                  data = datamc[, 2:nCol], outlvarz, az, wr, increment = 0)
    for (j in as.numeric(matz$.id)) {
      if (sum(az[[as.character(j)]] == 0) > 0) {
        datamc_imputed[j, 2:nCol] <- as.matrix(matz[matz$.id == j, -which(colnames(matz) == ".id")])
      }
    }
  }


  ## Update residuals and derive case weights with Hample weight function
  fitted_values <- datamc_imputed[, -1] %*% matrix(betas, ncol = 1) + intercept
  r <- datamc_imputed[, 1] - fitted_values
  r <- scaleResidualsByMAD(r)
  wr <- HampelWeightFunction(r, q1 = qnorm(0.95), q2 = qnorm(0.975), q3 = qnorm(0.999))
  names(wr) <- 1:nRow


  ## Loop until vector of regression coefficients has converged
  betas_old <- betas
  loop_counter <- 1
  difference <- 1
  difference_vec <- c()
  early_stop <- FALSE

  w0 <- which(wr < 1e-6)
  if (length(w0) > 0) {
    modelindex <- (1:nRow)[-w0]
  } else {
    modelindex <- 1:nRow
  }

  while (difference > tolerance & loop_counter <= maxiter & !early_stop) {
    if (verbose) {cat(paste0(" - Loop ", loop_counter, "..."))}


    ## Estimate linear model on weighted and imputed data
    modeldata <- diag(sqrt(wr[modelindex])) %*% datamc_imputed[modelindex, ]
    res.lm <- lm(formula, data = as.data.frame(modeldata))
    intercept <- ifelse(intercept_flag, res.lm$coefficients[1], 0)
    if (intercept_flag) {
      betas <- coef(res.lm)[-1]
    } else {
      betas <- coef(res.lm)
    }


    ## Calculate residuals of the imputed (not weighted data) and derive case weights
    fitted_values <- datamc_imputed[, -1] %*% matrix(betas, ncol = 1) + intercept
    r <- datamc_imputed[, 1] - fitted_values
    r <- scaleResidualsByMAD(r)
    wr <- HampelWeightFunction(r, q1 = qnorm(0.95), q2 = qnorm(0.975), q3 = qnorm(0.999))
    names(wr) <- 1:nRow


    ## Find cells contributing to outlyingness
    if (any(wr < 1)) {
      outliers <- as.list(which(wr < 1)) # casewise outliers in X-space (i.e. cases with wr[i] < 1)
      suppressWarnings(
        outlvars <- lapply(outliers, function (x) spadimo(data = datamc_imputed[, 2:nCol],
                                                          weights = wr, obs = x,
                                                          control = list(scaleFun  = Qn,
                                                                         nlatent   = 1,
                                                                         etas      = spadieta,
                                                                         csqcritv  = 0.975,
                                                                         stopearly = FALSE,
                                                                         trace     = FALSE,
                                                                         plot      = FALSE)))
      )
      outlvarz <- lapply(outlvars, function (x) {x$outlvars})
      az <- lapply(outlvars,function (x) {x$a})
      ob <- lapply(outlvars,function (x) {x$o.before})
      oa <- lapply(outlvars,function (x) {x$o.after})
      outliersx <- as.list(as.numeric(names(which(unlist(ob) > unlist(oa) * outlyingness.factor))))
      names(outliersx) <- as.character(unlist(outliersx))


      ## Update imputed data matrix
      matz <- ldply(outliersx, impute_outlying_cells,
                    data = datamc_imputed[, 2:nCol], outlvarz, az, wr, increment = 0)
      for (j in as.numeric(matz$.id)) {
        if (sum(az[[as.character(j)]] == 0) > 0) {
          datamc_imputed[j, 2:nCol] <- as.matrix(matz[matz$.id == j, -which(colnames(matz) == ".id")])
        }
      }
    }


    ## Update residuals and derive case weights with Hample weight function
    fitted_values <- datamc_imputed[, -1] %*% matrix(betas, ncol = 1) + intercept
    r <- datamc_imputed[, 1] - fitted_values
    r <- scaleResidualsByMAD(r)
    wr <- HampelWeightFunction(r, q1 = qnorm(0.95), q2 = qnorm(0.975), q3 = qnorm(0.999))
    names(wr) <- 1:nRow

    w0 <- which(wr < 1e-6)
    if (length(w0) > 0) {
      modelindex <- (1:nRow)[-w0]
    } else {
      modelindex <- 1:nRow
    }

    difference <- mean(abs(betas - betas_old))
    difference_vec <- c(difference, difference_vec)
    betas_old <- betas

    if (verbose) {cat(paste0("\r - Loop ", loop_counter,
                             " (difference = ", round(difference, 8), ")\n"))}
    if (length(difference_vec) >= 5) {
      early_stop <- all(abs(diff(difference_vec[1:5])) < sqrt(.Machine$double.eps))
    }
    loop_counter <- loop_counter + 1
  } # end of while loop


  if (difference > tolerance) {
    warning(paste("Method did not converge.",
                  "The scaled difference between the coefficient vectors is",
                  round(difference, digits = 4)))
  }


  ## Scale coefficients back to original scale
  if (verbose) {cat(" - Scale coefficients back to original scale...\n")}
  slopes <- attr(datamc, "Scale")[1] / attr(datamc, "Scale")[2:nCol] * betas
  if (intercept_flag) {
    if (center == "mean") {
      intercept <- mean(datam[, 1] - datam[, 2:nCol] %*% slopes)
    } else {
      intercept <- median(datam[, 1] - datam[, 2:nCol] %*% slopes)
    }
    coefs <- c(intercept, slopes)
    names(coefs)[1] <- "(Intercept)"
  } else {
    coefs <- slopes
  }


  ## Derive cells contributing to outlyingness
  cellwiseoutliers <- as.matrix(datamc[, -1] - datamc_imputed[, -1])


  ## Derive casewise outliers
  casewiseoutliers <- apply(cellwiseoutliers, 1, function (x) any(x != 0))


  ## Collect results
  data_imputed <- scale(datamc_imputed,
                        center = -attr(datamc, "Center") / attr(datamc, "Scale"),
                        scale = 1 / attr(datamc, "Scale"))
  yfit <- as.numeric(datam[, 2:nCol] %*% slopes + intercept)
  resid <- as.numeric(datam[, 1] - yfit)

  names(yfit) <- rownames(datam)
  names(resid) <- rownames(datam)
  names(wr) <- rownames(datam)
  names(casewiseoutliers) <- rownames(datam)
  rownames(cellwiseoutliers) <- rownames(datam)
  colnames(cellwiseoutliers) <- colnames(datam)[-1]


  ## End timer
  t_end <- proc.time() - t_start


  ## Output
  output <- list(coefficients     = coefs,
                 fitted.values    = yfit,
                 residuals        = resid,
                 weights          = wr,
                 data.imputed     = as.data.frame(data_imputed),
                 casewiseoutliers = casewiseoutliers,
                 cellwiseoutliers = cellwiseoutliers,
                 terms            = mt,
                 call             = call,
                 inputs           = inputs,
                 numloops         = loop_counter - 1,
                 time             = round(t_end[3], 3))
  class(output) <- "crm"
  return(output)
}
