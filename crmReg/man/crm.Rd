\name{crm}

\alias{crm}

\title{Cellwise Robust M-regression}

\description{Cellwise Robust M-regression.}

\usage{crm(formula, data, maxiter = 100, tolerance = 0.01, outlyingness.factor = 1,
    spadieta = seq(0.9, 0.1, -0.1), center = "median", scale = "qn",
    regtype = "LTS", alphaLTS = 0.5, seed = NULL, verbose = TRUE)}

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
\item{regtype}{type of robust regression. Choices are \code{"LTS"} (default) or \code{"MM"}.}
\item{alphaLTS}{parameter used by LTS regression. The percentage (roughly) of squared residuals whose sum will be minimized (default is \code{0.5}).}
\item{seed}{initial seed for random generator, like \link{.Random.seed} (default is \code{NULL}).}
\item{verbose}{should output be shown during the process (default is \code{TRUE}).}
}

\details{...}

\value{
\code{crm} returns a list object of class "crm" containing the following elements:
\item{coefficients}{a named vector of fitted coefficients.}
\item{fitted.values}{the fitted response values.}
\item{residuals}{the residuals, that is response minus fitted values.}
\item{weights}{the (case) weights of the residuals.}
\item{data.imputed}{the data as imputed by CRM.}
\item{casewiseoutliers}{a vector that indicates the casewise outliers with \code{TRUE} or \code{FALSE}.}
\item{cellwiseoutliers}{a matrix that indicates the cellwise outliers with 0, 1 or -1,
      depending on whether the variables are outlying upwards or downwards.}
\item{terms}{the terms object used.}
\item{call}{the matched call.}
\item{inputs}{the list of supplied input arguments.}
\item{numloops}{the number of iterations.}
\item{time}{the number of seconds passed to execute the CRM algorithm.}
}

\references{...}

\author{Peter Filzmoser, Sebastiaan Höppner, Iréne Ortner, Sven Serneels and Tim Verdonck}


\seealso{\code{\link{daprpr}}, \code{\link{spadimo}}}

\examples{
# Provide an example (TODO).
}
