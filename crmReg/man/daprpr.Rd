\name{daprpr}

\alias{daprpr}

\title{Data Preprocessing}

\description{Data preprocessing, classical and robus scaling.}

\usage{daprpr(Data, center.type, scale.type)}

\arguments{
\item{Data}{the data.}
\item{center.type}{type of centering as R function name (e.g. \code{"mean"}, \code{"median"}, \code{"l1median"}).}
\item{scale.type}{type of scaling as R function name (e.g. \code{"sd"}, \code{"qn"}, \code{"Sn"}, \code{"scaleTau2"}).}
}

\details{Given \code{center.type = "mean"} and \code{scale.type = "sd"}, \code{daprpr(...)} is equivalent to \code{scale(Data, center = TRUE, scale = TRUE)}.}

\value{
\code{daprpr} returns the scaled data with attributes \code{"Center"}, \code{"Scale"} and \code{"Type"}.}

\author{Sven Serneels}

\examples{
# Provide an example (TODO).
}
