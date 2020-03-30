\name{spadimo}

\alias{spadimo}

\title{SPArse DIrections of Maximal Outlyingness}

\description{Computes the sparse directions of maximal outlyings of a given
observation and shows diagnostic plots for analyzing that observation.}

\usage{spadimo(data, weights, obs,
        control = list(scaleFun = Qn, nlatent = 1, etas = NULL, csqcritv  = 0.975,
                       stopearly = FALSE, trace = FALSE, plot = TRUE))}

\arguments{
\item{data}{the data as a data frame.}
\item{weights}{a numeric vector containing the case weights from a robust estimator.}
\item{obs}{the (integer) case number under consideration.}
\item{control}{a list of options that control details of the \code{crm} algorithm. The following options are available:
\itemize{
\item{\code{scaleFun} \cr
function used for robust scaling the variables (e.g. \code{Qn}, \code{mad}, etc.).}

\item{\code{nlatent} \cr
integer number of latent variables for sparse PLS regression (via SNIPLS) (default is \code{1}).}

\item{\code{etas} \cr
vector of decreasing sparsity parameters (default is \code{NULL} in which case \code{etas = seq(0.9, 0.1, -0.05)} if n > p, otherwise \code{etas = seq(0.6, 0.1, -0.05)}).}

\item{\code{csqcritv} \cr
probability level for internal chi-squared quantile (used when n > p) (default is \code{0.975}).}

\item{\code{stopearly} \cr
if \code{TRUE}, method stops as soon as the reduced case is no longer outlying, else if \code{FALSE} (default) it loops through all values of eta.}

\item{\code{trace} \cr
should intermediate results be printed (default is \code{FALSE}).}

\item{\code{plot} \cr
should heatmaps and graph of the results be shown (default is \code{TRUE}).}
}
}
}

\details{
Given an observation that has been detected as an outlier, SPADIMO (Debruyne et al., 2019) finds
the subset of variables contributing most the outlier’s outlyingness. Here, the outlyingness of a
data point is defined as its robust Mahalanobis distance. The relevant variables are found by
checking the direction in which the observation is most outlying. SPADIMO estimates this direction
of maximal outlyingness in a sparse manner. Thereby, the method helps to understand in which way
an outlier lies out.
}

\value{
\code{spadimo} returns a list containing the following elements:
\item{outlvars}{vector containing individual variable names contributing most to \code{obs}'s outlyingness.}
\item{outlvarslist}{list of variables contributing to \code{obs}'s outlyingness for different values of \code{eta}.}
\item{a}{vector, the sparse direction of maximal outlyingness.}
\item{alist}{list of sparse directions of maximal outlyingness for different values of \code{eta}.}
\item{o.before}{outlyingness of original case (n < p) or PCA outlier flag (n >= p) before removing outlying variables.}
\item{o.after}{outlyingness of reduced case (n > p) or PCA outlier flag (n >= p) after removing outlying variables.}
\item{eta}{cutoff where \code{obs} is no longer outlying.}
\item{time}{time to execute the SPADIMO algorithm.}
\item{control}{a list with control parameters that are used.}
}

\references{
Debruyne, M., Hoppner, S., Serneels, S., and Verdonck, T. (2019). Outlyingness: Which variables contribute most?
\emph{Statistics and Computing}, 29 (4), 707–723.
\emph{DOI:10.1007/s11222-018-9831-5}
}

\author{Michiel Debruyne, Sebastiaan Hoppner, Sven Serneels, and Tim Verdonck}

\seealso{\code{\link{crm}}, \code{\link{predict.crm}}, \code{\link{cellwiseheatmap}}, \code{\link{daprpr}}}

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
}
