\name{cellwiseheatmap}

\alias{cellwiseheatmap}

\title{Heatmap of cellwise outliers}

\description{Heatmap of cellwise outliers.}

\usage{
cellwiseheatmap(cellwiseoutliers, data,
                col = c("blue", "lightgray", "red"),
                notecol.outlier = "white", notecol.clean = "black", notecex = 1,
                margins = c(9.5, 14), lhei = c(0.5, 15), lwid = c(0.1, 3.5),
                sepcolor = "white", sepwidth = c(0.01, 0.01))
}

\arguments{
\item{cellwiseoutliers}{a matrix that indicates the cellwise outliers with 0, 1 or -1,
      depending on whether the variables are outlying upwards or downwards.}
\item{data}{the data as a data frame that is shown in the cells, including row and column names.}
\item{col}{vector of colors used for downward outliers, clean cells and upward outliers respectively (default is \code{c("blue", "lightgray", "red")}).}
\item{notecol.outlier}{character string specifying the color for cellnote text of cellwise outliers (default is \code{"white"}).}
\item{notecol.clean}{character string specifying the color for cellnote text of clean cells (default is \code{"black"}).}
\item{notecex}{numeric scaling factor for cellnotes (default is \code{1}).}
\item{margins}{numeric vector of length 2 containing the margins (see \code{\link{par}(mar= *)}) for column and row names, respectively (default is \code{c(9.5, 14)}).}
\item{lhei}{numeric vector of length 2 containing the row height (default is \code{c(1, 15)}).}
\item{lwid}{numeric vector of length 2 containing the row width (default is \code{c(0.7, 3.5)}).}
\item{sepcolor}{character string specifying the color between the cells (default is \code{"white"}).}
\item{sepwidth}{vector of length 2 giving the width and height of the separator box drawn between the cells (default is \code{c(0.01, 0.01)}).}
}

\details{\code{cellwiseheatmap} plots a heatmap of cellwise outliers.}

\references{...}

\author{Peter Filzmoser, Sebastiaan Höppner, Iréne Ortner, Sven Serneels and Tim Verdonck}

\seealso{\code{\link{crm}}}

\examples{
# Provide an example (TODO).
}
