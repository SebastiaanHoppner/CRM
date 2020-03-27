\name{topgear}
\docType{data}
\alias{topgear}
\title{Top Gear car data}
\description{
The data set contains information on cars featured on the website of the popular BBC television
show \emph{Top Gear}. The original, full data set is available in the package \code{robustHD}.
}
\usage{data(topgear)}
\format{
A data frame containing 245 observations and 11 variables.
  \describe{
    \item{\code{log(Price)}}{the natural logarithm of the list price (in UK pounds)}
    \item{\code{log(Displacement)}}{the natural logarithm of the displacement of the engine (in cc).}
    \item{\code{log(BHP)}}{the natural logrithm of the power of the engine (in bhp).}
    \item{\code{log(Torque)}}{the natural logarithm of the torque of the engine (in lb/ft).}
    \item{\code{Acceleration}}{the time it takes the car to get from 0 to 62 mph (in seconds).}
    \item{\code{log(TopSpeed)}}{the natural logarithm of the car's top speed (in mph).}
    \item{\code{MPG}}{the combined fuel consuption (urban + extra urban; in miles per gallon).}
    \item{\code{Weight}}{the car's curb weight (in kg).}
    \item{\code{Length}}{the car's length (in mm).}
    \item{\code{Width}}{the car's width (in mm).}
    \item{\code{Height}}{the car's height (in mm).}
  }


}
\source{The original data set is available in the package \code{robustHD}.
The data were scraped from http://www.topgear.com/uk/ on 2014-02-24.}
\keyword{datasets}
\examples{
data(topgear)
str(topgear)
head(topgear)
summary(topgear)
}
