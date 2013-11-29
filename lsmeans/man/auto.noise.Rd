\name{auto.noise}
\alias{auto.noise}
\docType{data}
\title{Auto Pollution Filter Noise}
\description{Three-factor experiment comparing pollution-filter noise for two filters, three sizes of cars, and two sides of the car.}
\usage{auto.noise}

\format{
  A data frame with 36 observations on the following 4 variables.
  \describe{
    \item{\code{noise}}{Noise level in decibels - a numeric vector.}
    \item{\code{size}}{The size of the vehicle - an ordered factor with levels \code{S}, \code{M}, \code{L}.}
    \item{\code{type}}{Type of anti-pollution filter - a factor with levels \code{Std} and \code{Octel}}
    \item{\code{side}}{The side of the car where measurement was taken -- a factor with levels \code{R} and \code{L}.}
  }
}

\details{
The data are from a statement by Texaco, Inc. to the Air and Water Pollution
Subcommittee of the Senate Public Works Committee on June 26, 1973.    
Mr. John McKinley, President of Texaco, cited an automobile filter developed
by Associated Octel Company as effective in reducing pollution. However, 
questions had been raised about the effects of filters on vehicle performance, 
fuel consumption, exhaust gas back pressure, and silencing. On the last 
question, he referred to the data included here as evidence that the silencing
properties of the Octel filter were at least equal to those of standard silencers.}

\source{Dataset was imported from the Data and Story Library - \url{http://lib.stat.cmu.edu/DASL/Datafiles/airpullutionfiltersdat.html} (sic). However,
the factor levels were assigned meaningful names, and the observations were sorted
in random order as if this were the run order of the experiment.
}

\references{
A.Y. Lewin and M.F. Shakun (1976) \emph{Policy Sciences: Methodology and Cases}. Pergammon Press. p.313.
}
\examples{
require(lsmeans)
noise.lm = lm(noise ~ size * type * side, data = auto.noise)
lsmip(noise.lm, type ~ size | side)
}
\keyword{datasets}
