# Cellwise Robust M-Regression (CRM) - Simulation Study -------------------------------------------



rm(list = ls())
# Setup -------------------------------------------------------------------------------------------

# simulation setting - - - - - - - - - - - - - - -
n <- 400               # number of cases
p <- 50                # number of predictor variables
pct_case_out <- 0.05   # percentage of casewise outliers
pct_cell_out <- 0.10   # percentage of cellwise outliers for each casewise outlier
nsims <- 100           # number of simulations


# CRM input parameters - - - - - - - - - - - - - -
maxiter             <- 100
tolerance           <- 0.01
outlyingness.factor <- 1.5
spadieta            <- seq(0.9, 0.1, -0.1)



# Load packages -----------------------------------------------------------------------------------
library(crm)
library(MASS)
library(cellWise)
library(robustbase)



# Start simulation procedure ----------------------------------------------------------------------
results_MAE       <- data.frame(matrix(nrow = nsims, ncol = 5))
results_MSEP      <- data.frame(matrix(nrow = nsims, ncol = 5))
results_RMSEI     <- data.frame(matrix(nrow = nsims, ncol = 2))
results_precision <- data.frame(Precision = rep(NA, nsims))
results_recall    <- data.frame(Recall = rep(NA, nsims))
results_time      <- data.frame(Time = rep(NA, nsims))

names(results_MAE)   <- c("CRM", "MM", "DDC-MM", "OLS", "DDC-OLS")
names(results_MSEP)  <- c("CRM", "MM", "DDC-MM", "OLS", "DDC-OLS")
names(results_RMSEI) <- c("CRM", "DDC")


t_start <- proc.time()

set.seed(2019)
for (k_sim in 1:nsims) {
  
  t_start_round <- proc.time()
  cat(paste0("\nSimulation ", k_sim, "/", nsims, " ----------------------------------\n"))
  
  
  
  # Generate clean sample -------------------------------------------------------------------------
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
  
  
  
  # Add contamination in design matrix ------------------------------------------------------------
  Xc <- X
  contamination <- colMeans(X) + 6 * apply(X, 2, sd)
  
  case_outliers <- sample(n, size = n * pct_case_out)
  outliers_mat_flag <- matrix(FALSE, nrow = n, ncol = p)
  
  for (i in case_outliers) {
    cell_outliers <- sample(p, size = p * pct_cell_out)
    Xc[i, cell_outliers] <- contamination[cell_outliers] + rnorm(length(cell_outliers))
    outliers_mat_flag[i, cell_outliers] <- TRUE
  }
  outliers_mat_index <- which(outliers_mat_flag, arr.ind = TRUE)
  
  
  
  # Apply DDC to contaminated predictor data ------------------------------------------------------
  cat(" DDC:")
  DDC_Ximp <- DDC(Xc)$Ximp
  
  
  
  # Collect data samples --------------------------------------------------------------------------
  data_clean        <- cbind.data.frame(y, X)
  data_contaminated <- cbind.data.frame(y, Xc)
  data_DDC          <- cbind.data.frame(y, DDC_Ximp)
  names(data_clean) <- names(data_contaminated) <- names(data_DDC) <- c("y", paste0("X", 1:p))
  
  
  
  # Fit regression models -------------------------------------------------------------------------
  crm_fit <- crm(formula   = y ~ .,
                 data      = data_contaminated,
                 maxiter   = maxiter,
                 tolerance = tolerance,
                 outlyingness.factor = outlyingness.factor,
                 spadieta  = spadieta,
                 center    = "median",
                 scale     = "qn",
                 regtype   = "MM",
                 verbose   = FALSE)
  
  crm_lts_fit <- crm(formula   = y ~ .,
                     data      = data_contaminated,
                     maxiter   = maxiter,
                     tolerance = tolerance,
                     outlyingness.factor = outlyingness.factor,
                     spadieta  = spadieta,
                     center    = "median",
                     scale     = "qn",
                     regtype   = "LTS",
                     alphaLTS  = 0.5,
                     verbose   = FALSE)
  
  mm_fit      <- lmrob(formula = y ~ ., data = data_contaminated, method = "MM")
  ddc_mm_fit  <- lmrob(formula = y ~ ., data = data_DDC, method = "MM")
  
  lts_fit     <- ltsReg(formula = y ~ ., data = data_contaminated, alpha = 0.5)
  ddc_lts_fit <- ltsReg(formula = y ~ ., data = data_DDC, alpha = 0.5)
  
  ols_fit     <- lm(formula = y ~ ., data = data_contaminated)
  ddc_ols_fit <- lm(formula = y ~ ., data = data_DDC)
  
  
  
  # Evaluate performance --------------------------------------------------------------------------
  
  # Mean Absolute Error - - - - - - - - - - - - -
  results_MAE$CRM[k_sim]       <- mean(abs(    crm_fit$coefficients - betas))
  results_MAE$MM[k_sim]        <- mean(abs(     mm_fit$coefficients - betas))
  results_MAE$`DDC-MM`[k_sim]  <- mean(abs( ddc_mm_fit$coefficients - betas))
  results_MAE$OLS[k_sim]       <- mean(abs(    ols_fit$coefficients - betas))
  results_MAE$`DDC-OLS`[k_sim] <- mean(abs(ddc_ols_fit$coefficients - betas))
  
  
  # Mean Squared Error of Prediction - - - - - - -
  results_MSEP$CRM[k_sim]       <- mean((    crm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
  results_MSEP$MM[k_sim]        <- mean((     mm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
  results_MSEP$`DDC-MM`[k_sim]  <- mean(( ddc_mm_fit$fitted.values - data_clean$y)[-case_outliers]^2)
  results_MSEP$OLS[k_sim]       <- mean((    ols_fit$fitted.values - data_clean$y)[-case_outliers]^2)
  results_MSEP$`DDC-OLS`[k_sim] <- mean((ddc_ols_fit$fitted.values - data_clean$y)[-case_outliers]^2)
  
  
  # Root Mean Squared Error of Imputation - - - -
  results_RMSEI$CRM[k_sim] <- sqrt(mean((as.matrix(data_clean[, -1] - crm_fit$data.imputed[, -1]))^2))
  results_RMSEI$DDC[k_sim] <- sqrt(mean((as.matrix(data_clean[, -1] - data_DDC[, -1]))^2))
  
  
  # Precision & recall of detected cells - - - - -
  CM <- table(Prediction = ifelse(c(crm_fit$cellwiseoutliers) != 0, TRUE, FALSE),
              Reference  = c(outliers_mat_flag))
  
  results_precision$Precision[k_sim] <- CM[2, 2] / (CM[2, 1] + CM[2, 2])
  results_recall$Recall[k_sim]       <- CM[2, 2] / (CM[1, 2] + CM[2, 2])
  
  
  # Execution time - - - - - - - - - - - - - - - -
  results_time$Time[k_sim] <- crm_fit$time
  
  
  t_end_round <- proc.time() - t_start_round
  cat(paste0("\n ", round(100 * k_sim / nsims), "% - ", round(t_end_round[3]), " seconds \n"))
}
t_end <- proc.time() - t_start
cat(paste("\nTime elapsed:", floor(t_end[3]/60), "minutes", round(t_end[3] %% 60), "seconds\n\n"))



# Study results -----------------------------------------------------------------------------------

# Mean Absolute Error - - - - - - - - - - - - - -
boxplot(results_MAE, ylab = "MAE", col = "lightblue",
        ylim = c(0, max(results_MAE)), las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(results_MAE), col = "red", pch = 18, cex = 1.7)
text(rep(0, 8), labels = round(colMeans(results_MAE), 4), cex = 1.4,
     font = ifelse(1:5 == which.min(colMeans(results_MAE)), 2, 1))

# Mean Squared Error of Prediction - - - - - - - -
boxplot(results_MSEP, ylab = "MSEP", col = "lightblue",
        ylim = c(-0.5, max(results_MSEP)), las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(results_MSEP), col = "red", pch = 18, cex = 1.7)
text(rep(-0.5, 8), labels = round(colMeans(results_MSEP), 4), cex = 1.4,
     font = ifelse(1:5 == which.min(colMeans(results_MSEP)), 2, 1))

# Root Mean Squared Error of Imputation - - - - -
par(mfrow = c(1, 2))
boxplot(results_RMSEI, ylab = "RMSEI", col = "lightblue",
        ylim = c(0, max(results_RMSEI)), las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(results_RMSEI), col = "red", pch = 18, cex = 1.7)
text(rep(0, 2), labels = round(colMeans(results_RMSEI), 4), cex = 1.4,
     font = ifelse(1:2 == which.min(colMeans(results_RMSEI)), 2, 1))

# Precision & Recall - - - - - - - - - - - - - - -
boxplot(cbind(results_precision, results_recall), col = "lightblue",
        ylim = c(0, 1), ylab = "", las = 1, cex.axis = 1.2, cex.lab = 1.3)
points(colMeans(cbind(results_precision, results_recall)), col = "red", pch = 18, cex = 1.7)
text(rep(0, 2), labels = round(colMeans(cbind(results_precision, results_recall)), 4), cex = 1.4)
par(mfrow = c(1, 1))

# Execution time - - - - - - - - - - - - - - - - -
boxplot(results_time, col = "lightblue", ylab = "time (sec.)", las = 1,
        ylim = c(min(results_time) - 1, max(results_time) + 1))
points(colMeans(results_time), col = "red", pch = 18, cex = 1.7)
text(min(results_time) - 1, labels = round(colMeans(results_time), 4), cex = 1.4)


