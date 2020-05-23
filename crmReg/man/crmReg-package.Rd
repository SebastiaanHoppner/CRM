\name{crmReg-package}
\alias{crmReg-package}
\alias{crmReg}
\docType{package}
\title{
Cellwise Robust M-regression and SPADIMO
}
\description{
Method for fitting a cellwise robust linear M-regression model (CRM, Filzmoser et al. (2020) <DOI:10.1016/j.csda.2020.106944>) that yields both
a map of cellwise outliers consistent with the linear model, and a vector of regression coefficients
that is robust against vertical outliers and leverage points. As a by-product, the method yields an
imputed data set that contains estimates of what the values in cellwise outliers would need to
amount to if they had fit the model. The package also provides diagnostic tools for analyzing
casewise and cellwise outliers using sparse directions of maximal outlyingness (SPADIMO, Debruyne et al. (2019) <DOI:10.1007/s11222-018-9831-5>).
}
\details{
\tabular{ll}{
Package: \tab crmReg \cr
Type: \tab Package \cr
Version: \tab 1.0.1 \cr
Date: \tab 2020-03-26 \cr
License: \tab GPL (>=2) \cr
}

The crmReg package provides the implementation of the Cellwise Robust M-regression (CRM) algorithm
(Filzmoser et al., 2020) and the SPArse DIrections of Maximal Outlyingness (SPADIMO) algorithm
(Debruyne et al., 2019). The package also includes a predict function for fitted CRM regression
models, a function for creating heatmaps of cellwise outliers, and a data preprocessing function
for centering and scaling the data as used by CRM.

Given an observation that has been detected as an outlier, SPADIMO (Debruyne et al., 2019) finds
the subset of variables contributing most the outlier’s outlyingness. Here, the outlyingness of a
data point is defined as its robust Mahalanobis distance. The relevant variables are found by
checking the direction in which the observation is most outlying. SPADIMO estimates this direction
of maximal outlyingness in a sparse manner. Thereby, the method helps to understand in which way
an outlier lies out.

The SPADIMO algorithm allows us to introduce the cellwise robust M-regression (CRM) estimator
(Filzmoser et al., 2020) as a linear regression estimator that intrinsically yields both a map of
cellwise outliers consistent with the linear model, and a vector of regression coefficients that
is robust against vertical outliers and leverage points.
As a by-product, the method yields a weighted and imputed data set that
contains estimates of what the values in cellwise outliers would need to amount to if they had fit
the model. The CRM method consists of an iteratively reweighted least squares procedure where
SPADIMO is applied at each iteration to detect the cells that contribute most to outlyingness.
As such, CRM detects deviating data cells consistent with a linear model.

The package contains five main functions.

The function \code{\link{spadimo}} computes the sparse directions of maximal outlyings of a given
observation and shows diagnostic plots for analyzing that observation.

The function \code{\link{crm}} fits a cellwise robust M-regression estimator. Besides a vector of
regression coefficients, the function returns an imputed data set that contains estimates of what
the values in cellwise outliers would need to amount to if they had fit the model. The output of
\code{crm} is a list object of class "\code{crm}".

The function \code{\link{predict.crm}} obtains predictions from a fitted object of class "\code{crm}".

The function \code{\link{cellwiseheatmap}} makes a heatmap of cellwise outliers which are typically
the result of a call to the \code{crm} function.

The function \code{\link{daprpr}} preprocesses the data by classical or robust centering and scaling.
}
\author{
Peter Filzmoser, Sebastiaan Hoppner, Irene Ortner, Sven Serneels, and Tim Verdonck

Maintainer: Sebastiaan Hoppner <sebastiaan.hoppner@kuleuven.be>
}
\references{
Debruyne, M., Hoppner, S., Serneels, S., and Verdonck, T. (2019). Outlyingness: Which variables contribute most?
\emph{Statistics and Computing}, 29 (4), 707–723.
\emph{DOI:10.1007/s11222-018-9831-5}

Filzmoser, P., Hoppner, S., Ortner, I., Serneels, S., and Verdonck, T. (2020). Cellwise Robust M regression.
\emph{Computational Statistics and Data Analysis}, 147, 106944.
\emph{DOI:10.1016/j.csda.2020.106944}
}

\seealso{
\code{\link{crm}}, \code{\link{spadimo}}, \code{\link{predict.crm}}, \code{\link{cellwiseheatmap}}, \code{\link{daprpr}}
}
\examples{
library(crmReg)
data(topgear)

# get case weights from a robust estimator (covMCD function in robustbase package):
MCD <- robustbase::covMcd(topgear, alpha = 0.5)

# SPADIMO with diagnostic plots:
# Example 1:
Peugeot <- spadimo(data = topgear,
                   weights = MCD$mcd.wt,
                   obs = which(rownames(topgear) == "Peugeot 107"))
# check the plots!
# individual variable names contributing most to Peugeot 107's outlyingness:
print(Peugeot$outlvars)
# sparse direction of maximal outlyingness with eta = Peugeot$eta:
print(Peugeot$a)
# default SPADIMO control parameters:
print(Peugeot$control)

# Example 2:
Bugatti <- spadimo(data = topgear,
                   weights = MCD$mcd.wt,
                   obs = which(rownames(topgear) == "Bugatti Veyron"),
                   control = list(stopearly = TRUE, trace = TRUE, plot = TRUE))
# check the plots!
# individual variable names contributing most to Bugatti Veyron's outlyingness:
print(Bugatti$outlvars)
# sparse direction of maximal outlyingness with eta = Bugatti$eta:
print(Bugatti$a)

# fit Cellwise Robust M-regression:
crmfit <- crm(formula = MPG ~ ., data = topgear)

# estimated regression coefficients and detected casewise outliers:
print(crmfit$coefficients)
print(rownames(topgear)[which(crmfit$casewiseoutliers)])

# fitted response values (MPG) versus true response values:
plot(topgear$MPG, crmfit$fitted.values, xlab = "True MPG", ylab = "Fitted MPG")
abline(a = 0, b = 1)

# residuals:
plot(crmfit$residuals, ylab = "Residuals")
text(x = which(crmfit$residuals > 30), y = crmfit$residuals[which(crmfit$residuals > 30)],
     labels = rownames(topgear)[which(crmfit$residuals > 30)], pos = 2)

print(cbind.data.frame(car = rownames(topgear),
                       MPG = topgear$MPG)[which(crmfit$residuals > 30), ])

# cellwise heatmap of casewise outliers:
cellwiseheatmap(cellwiseoutliers = crmfit$cellwiseoutliers[which(crmfit$casewiseoutliers), ],
                data = round(topgear[which(crmfit$casewiseoutliers), -7], 2),
                col.scale.factor = 1/4)
# check the plotted heatmap!
}
