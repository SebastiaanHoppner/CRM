# Nutrient contents -------------------------------------------------------------------------------



# Load & preprocess data --------------------------------------------------------------------------
rm(list = ls())
data(nutrients, package = "robCompositions")
nutrients <- as.data.frame(nutrients)

# change rownames to the English name of the food product combined with the row numbers:
rownames(nutrients) <- paste(nutrients$name_E, 1:nrow(nutrients))

# select numerical variables:
nutrients <- nutrients[1:200, c(23, 14, 15, 17, 18, 20)]

# check for missing values:
print(sum(rowSums(is.na(nutrients)) > 0))
print(apply(nutrients, 2, function (x) sum(is.na(x))))
# 7 out of 200 food products have at least 1 missing value

nutrients <- na.omit(nutrients) # remove food products with at least 1 NA value
print(dim(nutrients))

# log-transform all variables:
transNutrients <- log(nutrients + 0.01)
colnames(transNutrients) <- paste0("log.", names(nutrients))
# transNutrients is the data set that we analyse with CRM.



# Cellwise Robust M-regression --------------------------------------------------------------------
library(crm)
help(crm)

crm_fit <- crm(formula   = log.cholesterol ~ .,
               data      = transNutrients,
               maxiter   = 100,
               tolerance = 0.01,
               outlyingness.factor = 1,
               spadieta  = seq(0.9, 0.1, -0.1),
               center    = "median",
               scale     = "sd",
               regtype   = "MM",
               seed      = 2019,
               verbose   = FALSE)



# Analyse results ---------------------------------------------------------------------------------

# execution time - - - - - - - - - - - - - - - - -
print(crm_fit$time)


# coefficients - - - - - - - - - - - - - - - - - -
print(round(crm_fit$coefficients, 5))


# casewise outliers - - - - - - - - - - - - - - -
casesOfInterest <- which(crm_fit$casewiseoutliers)
print(length(casesOfInterest))
print(rownames(nutrients)[casesOfInterest])
# 26 out of 193 food products are considered as casewise outliers by the CRM algorithm


# heatmap of cellwise outliers - - - - - - - - - -

# original nutrients data:
nutrients_casesOfInterest <- round(nutrients[casesOfInterest, -1], 1)
rownames(nutrients_casesOfInterest) <- substr(rownames(nutrients_casesOfInterest), 1, 
                                              nchar(rownames(nutrients_casesOfInterest)) - c(2, rep(3, length(casesOfInterest)-1)))

cellwiseheatmap(cellwiseoutliers = crm_fit$cellwiseoutliers[casesOfInterest, ],
                data = nutrients_casesOfInterest,
                margins = c(10, 22), notecex = 2)

# imputed nutrients data by CRM:
nutrients.imputed <- exp(crm_fit$data.imputed)
colnames(nutrients.imputed) <- colnames(nutrients)
nutrients.imputed_casesOfInterest <- round(nutrients.imputed[casesOfInterest, -1], 1)
rownames(nutrients.imputed_casesOfInterest) <- rownames(nutrients_casesOfInterest)

cellwiseheatmap(cellwiseoutliers = crm_fit$cellwiseoutliers[casesOfInterest, ],
                data = nutrients.imputed_casesOfInterest,
                margins = c(10, 22), notecex = 2)



# Cross-validation --------------------------------------------------------------------------------
library(caret)
library(robustbase)

nfolds <- 10
set.seed(2019)
folds <- createFolds(y = transNutrients$log.cholesterol, k = nfolds)

RMSE <- data.frame(CRM = rep(NA, nfolds), MM = rep(NA, nfolds))

set.seed(2019)
for (k in 1:nfolds) {
  cat(paste("\n Fold", k, "out of", nfolds, "\n"))
  
  train <- transNutrients[-folds[[k]], ]
  test  <- transNutrients[ folds[[k]], ]
  
  crm_fit <- crm(formula   = log.cholesterol ~ .,
                 data      = train,
                 maxiter   = 100,
                 tolerance = 0.01,
                 outlyingness.factor = 1,
                 spadieta  = seq(0.9, 0.1, -0.1),
                 center    = "median",
                 scale     = "qn",
                 regtype   = "MM",
                 seed      = 2019,
                 verbose   = FALSE)
  
  mm_fit <- lmrob(formula = log.cholesterol ~ ., data = train, method = "MM")
  
  pred_test_crm <- predict(crm_fit, newdata = test)
  pred_test_mm  <- predict(mm_fit,  newdata = test)
  
  RMSE$CRM[k] <- sqrt(mean((test$log.cholesterol - pred_test_crm)^2, trim = 0.1))
  RMSE$MM[k]  <- sqrt(mean((test$log.cholesterol - pred_test_mm)^2,  trim = 0.1))
}

boxplot(RMSE, ylab = "10% trimmed RMSEP", col = "olivedrab4",
        ylim = c(0, max(RMSE)), cex.lab = 1.4, cex.axis = 1.4, cex.main = 2)
points(colMeans(RMSE), pch = 18, cex = 1.5)
text(rep(0, 2), labels = round(colMeans(RMSE), 4), cex = 2,
     font = ifelse(1:2 == which.min(colMeans(RMSE)), 2, 1))


