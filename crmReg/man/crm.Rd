\name{crm}

\alias{crm}

\title{Cellwise Robust M-regression}

\description{Fits a cellwise robust M-regression estimator. Besides a vector of
regression coefficients, the function returns an imputed data set that contains estimates of what
the values in cellwise outliers would need to amount to if they had fit the model.}

\usage{crm(formula, data, maxiter = 100, tolerance = 0.01, outlyingness.factor = 1,
    spadieta = seq(0.9, 0.1, -0.1), center = "median", scale = "qn",
    regtype = "MM", alphaLTS = NULL, seed = NULL, verbose = TRUE)}

\arguments{
\item{formula}{an lm-style formula object specifying which relationship to estimate.}
\item{data}{the data as a data frame.}
\item{maxiter}{maximum number of iterations (default is \code{100}).}
\item{tolerance}{obtain optimal regression coefficients to within a certain tolerance (default is \code{0.01}).}
\item{outlyingness.factor}{numeric value, larger or equal to 1 (default). Only cells are altered of
cases for which the original outlyingness (before SPADIMO) is larger than
outlyingness.factor * outlyingness AFTER SPADIMO. The larger this factor, the fewer cells are imputed.}
\item{spadieta}{the sparsity parameter to start internal outlying cell detection with, must be in the range [0,1] (default is \code{seq(0.9, 0.1, -0.1)}).}
\item{center}{how to center the data. A string that matches the R function to be used for centering (default is \code{"median"}).}
\item{scale}{how to scale the data. Choices are "no" (no scaling) or a string matching the R function to be used for scaling (default is \code{"qn"}).}
\item{regtype}{type of robust regression. Choices are \code{"MM"} (default) or \code{"LTS"}.}
\item{alphaLTS}{parameter used by LTS regression. The percentage (roughly) of squared residuals whose sum will be minimized (default is \code{0.5}).}
\item{seed}{initial seed for random generator, like \link{.Random.seed} (default is \code{NULL}).}
\item{verbose}{should output be shown during the process (default is \code{TRUE}).}
}

\details{The cellwise robust M-regression (CRM) estimator (Filzmoser et al., 2020) is a linear
regression estimator that intrinsically yields both a map of cellwise outliers consistent with the
linear model, and a vector of regression coefficients that is robust against vertical outliers and leverage points.
As a by-product, the method yields a weighted and imputed data set that contains estimates of what
the values in cellwise outliers would need to amount to if they had fit the model. The CRM method
consists of an iteratively reweighted least squares procedure where SPADIMO is applied at each
iteration to detect the cells that contribute most to outlyingness. As such, CRM detects deviating
data cells consistent with a linear model.}

\value{
\code{crm} returns a list object of class "\code{crm}" containing the following elements:
\item{coefficients}{a named vector of fitted coefficients.}
\item{fitted.values}{the fitted response values.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{weights}{the (case) weights of the residuals.}
\item{data.imputed}{the data as imputed by CRM.}
\item{casewiseoutliers}{a vector that indicates the casewise outliers with \code{TRUE} or \code{FALSE}.}
\item{cellwiseoutliers}{a matrix that indicates the cellwise outliers as the (scaled) difference between the original data and imputed data, both scaled and centered.}
\item{terms}{the terms object used.}
\item{call}{the matched call.}
\item{inputs}{the list of supplied input arguments.}
\item{numloops}{the number of iterations.}
\item{time}{the number of seconds passed to execute the CRM algorithm.}
}

\references{Filzmoser, P., Hoeppner, S., Ortner, I., Serneels, S., and Verdonck, T. (2020). Cellwise Robust M regression.
\emph{Computational Statistics and Data Analysis}, 147(106944). DOI: 10.1016/j.csda.2020.106944}

\author{Peter Filzmoser, Sebastiaan Hoeppner, Iréne Ortner, Sven Serneels, and Tim Verdonck}


\seealso{\code{\link{spadimo}}, \code{\link{predict.crm}}, \code{\link{cellwiseheatmap}}, \code{\link{daprpr}}}

\examples{
library(crmReg)
data(topgear)

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
