\name{predict.crm}

\alias{predict.crm}

\title{Predict method for CRM fits}

\description{Obtains predictions from a fitted \code{crm} object.}

\usage{
\method{predict}{crm}(object, newdata = NULL, \dots)
}

\arguments{
\item{object}{a fitted object of class "\code{crm}".}
\item{newdata}{optionally, a data frame in which to look for variables with which to predict.
If omitted, the fitted coefficients are used.}
\item{...}{further arguments passed to or from other methods.}
}

\details{\code{predict.crm} produces predicted values, obtained by evaluating the fitted
\code{\link{crm}} object on the data frame \code{newdata}.}

\value{
\code{predict.crm} returns a vector of predicted response values.
}

\references{
Filzmoser, P., Hoppner, S., Ortner, I., Serneels, S., and Verdonck, T. (2020). Cellwise Robust M regression.
\emph{Computational Statistics and Data Analysis}.
\emph{DOI:10.1016/j.csda.2020.106944}
}

\author{Peter Filzmoser, Sebastiaan Hoppner, Irene Ortner, Sven Serneels, and Tim Verdonck}

\seealso{\code{\link{crm}}, \code{\link{spadimo}}, \code{\link{cellwiseheatmap}}, \code{\link{daprpr}}}

\examples{
library(crmReg)
data(topgear)

train <- topgear[1:200, ]
test <- topgear[201:245, ]

crmfit <- crm(formula = MPG ~ ., data = train, seed = 2020)

estimated_MPG_test <- predict(crmfit, newdata = test)

plot(test$MPG, estimated_MPG_test, xlab = "True MPG", ylab = "Estimated MPG")
abline(a = 0, b = 1)
}
