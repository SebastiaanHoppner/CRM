\name{crmReg-package}
\alias{crmReg-package}
\alias{crmReg}
\docType{package}
\title{
Cellwise Robust M-regression and SPADIMO
}
\description{
Cellwise Robust M-regression (CRM) and SPArse DIrections of Maximal Outlyiness (SPADIMO).
}
\details{
\tabular{ll}{
Package: \tab crmReg \cr
Type: \tab Package \cr
Version: \tab 1.0.0 \cr
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
Peter Filzmoser, Sebastiaan Höppner, Iréne Ortner, Sven Serneels, and Tim Verdonck

Maintainer: Sebastiaan Höppner <sebastiaan.hoppner@kuleuven.be>
}
\references{
Debruyne, M., Höppner, S., Serneels, S., and Verdonck, T. (2019). Outlyingness: Which variables contribute mots?
\emph{Statistics and Computing}, 29 (4), 707–723. DOI: 10.1007/s11222-018-9831-5

Filzmoser, P., Höppner, S., Ortner, I., Serneels, S., and Verdonck, T. (2020). Cellwise Robust M regression.
\emph{Computational Statistics and Data Analysis}. DOI: 10.1016/j.csda.2020.106944
}

\seealso{
\code{\link{crm}}, \code{\link{spadimo}}, \code{\link{predict.crm}}, \code{\link{cellwiseheatmap}}, \code{\link{daprpr}}
}
\examples{
# Provide an example (TODO).
}
