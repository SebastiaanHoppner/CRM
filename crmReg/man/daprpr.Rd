\name{daprpr}

\alias{daprpr}

\title{Data Preprocessing}

\description{Data preprocessing, classical and robust centering and scaling.}

\usage{daprpr(Data, center.type, scale.type)}

\arguments{
\item{Data}{the data.}
\item{center.type}{type of centering as R function name (e.g. \code{"mean"}, \code{"median"}, \code{"l1median"}).}
\item{scale.type}{type of scaling as R function name (e.g. \code{"sd"}, \code{"qn"}, \code{"Sn"}, \code{"scaleTau2"}).}
}

\details{\code{daprpr} preprocesses the data by classical or robust centering and scaling.
Given \code{center.type = "mean"} and \code{scale.type = "sd"}, function \code{daprpr} is
equivalent to \code{scale(Data, center = TRUE, scale = TRUE)}.}

\value{
\code{daprpr} returns the scaled data with attributes \code{"Center"}, \code{"Scale"} and \code{"Type"}.}

\author{Sven Serneels}

\seealso{\code{\link{crm}}, \code{\link{spadimo}}, \code{\link{predict.crm}}, \code{\link{cellwiseheatmap}}}

\examples{
library(crmReg)
data(topgear)

topgear_centered_scaled <- daprpr(topgear, center.type = "median", scale.type = "qn")

boxplot(topgear_centered_scaled)
attributes(topgear_centered_scaled)$Type
attributes(topgear_centered_scaled)$Center
attributes(topgear_centered_scaled)$Scale
}
