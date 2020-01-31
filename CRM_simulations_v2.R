# Cellwise Robust M-regression (CRM) - Simulation Study -------------------------------------------



rm(list = ls())
# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400               # number of cases
p <- 50                # number of predictor variables
pct_case_out <- 0.05   # percentage of casewise outliers
pct_cell_out <- 0.10   # percentage of cellwise outliers for each casewise outlier
nsims <- 10            # number of simulations for each value of k
k_seq <- 0:8           # sequence of parameter values that determine the contamination


# CRM input parameters - - - - - - - - - - - - - -
maxiter             <- 100
tolerance           <- 0.01
outlyingness.factor <- 1.5
spadieta            <- seq(0.9, 0.1, -0.1)



# Load packages -----------------------------------------------------------------------------------
library(crm)
library(MASS)
library(cellWise)
library(lubridate)
library(robustbase)



# Start simulation procedure ----------------------------------------------------------------------
results_MAE       <- list()
results_MSEP      <- list()
results_RMSEI     <- list()
results_precision <- list()
results_recall    <- list()
results_time      <- list()

t_start <- proc.time()
set.seed(2020)


for (k_value in k_seq) {
  
  k_index <- which(k_seq == k_value)
  cat(paste("\n*", k_index, "out of", length(k_seq), "========================================\n"))
  
  results_MAE_k       <- data.frame(matrix(nrow = nsims, ncol = 5))
  results_MSEP_k      <- data.frame(matrix(nrow = nsims, ncol = 5))
  results_RMSEI_k     <- data.frame(matrix(nrow = nsims, ncol = 2))
  results_precision_k <- data.frame(Precision = rep(NA, nsims))
  results_recall_k    <- data.frame(Recall = rep(NA, nsims))
  results_time_k      <- data.frame(Time = rep(NA, nsims))
  
  names(results_MAE_k)   <- c("CRM", "MM", "DDC-MM", "OLS", "DDC-OLS")
  names(results_MSEP_k)  <- c("CRM", "MM", "DDC-MM", "OLS", "DDC-OLS")
  names(results_RMSEI_k) <- c("CRM", "DDC")
  
  
  for (j_sim in 1:nsims) {
    
    
    # Generate clean sample -----------------------------------------------------------------------
    mu <- rep(0, p)
    Sigma <- diag(p)
    Sigma[(row(Sigma) - col(Sigma)) == 1 | row(Sigma) - col(Sigma) == -1] <- 0.5
    X <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
    
    slopes <- rnorm(n = p, mean = 0, sd = 1)
    slopes <- 10 * slopes / sqrt(sum(slopes^2))
    intercept <- 10
    noise <- rnorm(n = n, mean = 0, sd = 0.5)
    
    y <- intercept + X %*% slopes + noise
    betas <- c(intercept, slopes)
    
    
    
    # Add contamination in design matrix ----------------------------------------------------------
    Xc <- X
    contamination <- colMeans(X) + k_value * apply(X, 2, sd)
    
    case_outliers <- sample(n, size = n * pct_case_out)
    outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)
    
    for (i in case_outliers) {
      cell_outliers <- sample(p, size = p * pct_cell_out)
      Xc[i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
      outliers_mat_flag[i, cell_outliers] <- TRUE
    }
    outliers_mat_index <- which(outliers_mat_flag, arr.ind = TRUE)
    
    
    
    # Apply DDC to contaminated predictor data ----------------------------------------------------
    cat(" DDC:")
    DDC_Ximp <- DDC(Xc)$Ximp
    
    
    
    # Collect data samples ------------------------------------------------------------------------
    data_clean        <- cbind.data.frame(y, X)
    data_contaminated <- cbind.data.frame(y, Xc)
    data_DDC          <- cbind.data.frame(y, DDC_Ximp)
    names(data_clean) <- names(data_contaminated) <- names(data_DDC) <- c("y", paste0("X", 1:p))
    
    
    
    # Fit regression models -----------------------------------------------------------------------
    crm_fit <- suppressWarnings(crm(formula   = y ~ .,
                                    data      = data_contaminated,
                                    maxiter   = maxiter,
                                    tolerance = tolerance,
                                    outlyingness.factor = outlyingness.factor,
                                    spadieta  = spadieta,
                                    center    = "median",
                                    scale     = "qn",
                                    regtype   = "MM",
                                    verbose   = FALSE))
    
    mm_fit      <- suppressWarnings(lmrob(formula = y ~ ., data = data_contaminated, method = "MM"))
    ddc_mm_fit  <- suppressWarnings(lmrob(formula = y ~ ., data = data_DDC, method = "MM"))
    
    ols_fit     <- lm(formula = y ~ ., data = data_contaminated)
    ddc_ols_fit <- lm(formula = y ~ ., data = data_DDC)
    
    
    
    # Evaluate performance ------------------------------------------------------------------------
    
    # Mean Absolute Error - - - - - - - - - - - - 
    results_MAE_k$`CRM`[j_sim]     <- mean(abs(    crm_fit$coefficients - betas))
    results_MAE_k$`MM`[j_sim]      <- mean(abs(     mm_fit$coefficients - betas))
    results_MAE_k$`DDC-MM`[j_sim]  <- mean(abs( ddc_mm_fit$coefficients - betas))
    results_MAE_k$`OLS`[j_sim]     <- mean(abs(    ols_fit$coefficients - betas))
    results_MAE_k$`DDC-OLS`[j_sim] <- mean(abs(ddc_ols_fit$coefficients - betas))
    
    
    # Mean Squared Error of Prediction - - - - - -
    results_MSEP_k$`CRM`[j_sim]     <- mean((    crm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
    results_MSEP_k$`MM`[j_sim]      <- mean((     mm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
    results_MSEP_k$`DDC-MM`[j_sim]  <- mean(( ddc_mm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
    results_MSEP_k$`OLS`[j_sim]     <- mean((    ols_fit$fitted.values - data_clean$y)[-case_outliers]^2)
    results_MSEP_k$`DDC-OLS`[j_sim] <- mean((ddc_ols_fit$fitted.values - data_clean$y)[-case_outliers]^2)
    
    
    # Root Mean Squared Error of Imputation - - - 
    results_RMSEI_k$CRM[j_sim] <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit$data.imputed[, -1]))^2))
    results_RMSEI_k$DDC[j_sim] <- sqrt(mean((as.matrix(data_clean[, -1] - data_DDC[, -1]))^2))
    
    
    # Precision & recall of detected cells - - - -
    CM <- table(Prediction = ifelse(c(crm_fit$cellwiseoutliers) != 0, TRUE, FALSE),
                Reference  = c(outliers_mat_flag))
    
    results_precision_k$Precision[j_sim] <- CM[2, 2] / (CM[2, 1] + CM[2, 2])
    results_recall_k$Recall[j_sim]       <- CM[2, 2] / (CM[1, 2] + CM[2, 2])
    
    
    # Execution time - - - - - - - - - - - - - - -
    results_time_k$Time[j_sim] <- crm_fit$time
    
    
    # Elapsed time - - - - - - - - - - - - - - - -
    t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
    cat(paste0(" - ", round(100 * (j_sim + (k_index - 1) * nsims) / (nsims * length(k_seq)), 2),
               "% - ", sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))
    
    
  } # end of forloop "j_sim in 1:nsims"
  
  
  results_MAE[[k_index]]       <- results_MAE_k
  results_MSEP[[k_index]]      <- results_MSEP_k
  results_RMSEI[[k_index]]     <- results_RMSEI_k
  results_precision[[k_index]] <- results_precision_k
  results_recall[[k_index]]    <- results_recall_k
  results_time[[k_index]]      <- results_time_k
  
  
} # end of forloop "k_value in k_seq"


t_end <- seconds_to_period(round((proc.time() - t_start)[3]))
cat(paste("\nTime elapsed:",
          sprintf("%02d:%02d:%02d", t_end@hour, minute(t_end), second(t_end)), "\n\n"))



# Study results -----------------------------------------------------------------------------------
results_MAE       <- cbind.data.frame(do.call(rbind, results_MAE),       k = sort(rep(k_seq, nsims)))
results_MSEP      <- cbind.data.frame(do.call(rbind, results_MSEP),      k = sort(rep(k_seq, nsims)))
results_RMSEI     <- cbind.data.frame(do.call(rbind, results_RMSEI),     k = sort(rep(k_seq, nsims)))
results_precision <- cbind.data.frame(do.call(rbind, results_precision), k = sort(rep(k_seq, nsims)))
results_recall    <- cbind.data.frame(do.call(rbind, results_recall),    k = sort(rep(k_seq, nsims)))
results_time      <- cbind.data.frame(do.call(rbind, results_time),      k = sort(rep(k_seq, nsims)))

results_MAE_mean       <- aggregate(. ~ k, data = results_MAE,       FUN = mean)
results_MSEP_mean      <- aggregate(. ~ k, data = results_MSEP,      FUN = mean)
results_RMSEI_mean     <- aggregate(. ~ k, data = results_RMSEI,     FUN = mean)
results_precision_mean <- aggregate(. ~ k, data = results_precision, FUN = mean)
results_recall_mean    <- aggregate(. ~ k, data = results_recall,    FUN = mean)
results_time_mean      <- aggregate(. ~ k, data = results_time,      FUN = mean)


# Mean Absolute Error - - - - - - - - - - - - - - 
plot(results_MAE_mean$k, results_MAE_mean$`CRM`, type = "l", xlab = "k", ylab = "Average MAE",
     ylim = range(results_MAE_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(results_MAE_mean$k, results_MAE_mean$`MM`,      lwd = 2, lty = 2)
lines(results_MAE_mean$k, results_MAE_mean$`DDC-MM`,  lwd = 2, lty = 3)
lines(results_MAE_mean$k, results_MAE_mean$`OLS`,     lwd = 2, lty = 4)
lines(results_MAE_mean$k, results_MAE_mean$`DDC-OLS`, lwd = 2, lty = 5)
legend("topleft", legend = names(results_MAE_mean)[-1], lty = 1:5, lwd = 2)


# Mean Squared Error of Prediction - - - - - - - -
plot(results_MSEP_mean$k, results_MSEP_mean$`CRM`, type = "l", xlab = "k", ylab = "Average MSEP",
     ylim = range(results_MSEP_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(results_MSEP_mean$k, results_MSEP_mean$`MM`,      lwd = 2, lty = 2)
lines(results_MSEP_mean$k, results_MSEP_mean$`DDC-MM`,  lwd = 2, lty = 3)
lines(results_MSEP_mean$k, results_MSEP_mean$`OLS`,     lwd = 2, lty = 4)
lines(results_MSEP_mean$k, results_MSEP_mean$`DDC-OLS`, lwd = 2, lty = 5)
legend("left", legend = names(results_MSEP_mean)[-1], lty = 1:5, lwd = 2)


# Root Mean Squared Error of Imputation - - - - - 
plot(results_RMSEI_mean$k, results_RMSEI_mean$CRM, type = "l", xlab = "k", ylab = "Average RMSEI",
     ylim = range(results_RMSEI_mean[, -1]), las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
lines(results_RMSEI_mean$k, results_RMSEI_mean$DDC, lwd = 2, lty = 2)
legend("left", legend = names(results_RMSEI_mean)[-1], lty = 1:2, lwd = 2)


# Precision detected cells - - - - - - - - - - - -
par(mfrow = c(2, 1))
plot(results_precision_mean, type = "l", xlab = "k", ylab = "Average precision",
     las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, ylim = c(0, 0.7))


# Recall of detected cells - - - - - - - - - - - -
plot(results_recall_mean, type = "l", xlab = "k", ylab = "Average recall",
     las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2, ylim = c(0, 1))
par(mfrow = c(1, 1))


# Execution time - - - - - - - - - - - - - - - - -
plot(results_time_mean, type = "l", xlab = "k", ylab = "Average time (sec.)",
     las = 1, cex.axis = 1.2, cex.lab = 1.3, lwd = 2)
cat("CRM average execution time:", round(mean(results_time$Time), 1), "seconds")


