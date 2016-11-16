%    Copyright (c) 2012-2016 Russell V. Lenth                                %

\name{fiber}
\alias{fiber}
\docType{data}
\title{
Fiber data
}
\description{
Fiber data from Montgomery Design (8th ed.), p.656 (Table 15.10). Useful as a simple analysis-of-covariance example.
}
\usage{fiber}
\format{
  A data frame with 15 observations on the following 3 variables.
  \describe{
    \item{\code{machine}}{a factor with levels \code{A} \code{B} \code{C}. The primary factor of interest.}
    \item{\code{strength}}{a numeric vector. The response variable.}
    \item{\code{diameter}}{a numeric vector. A covariate.}
  }
}
\details{
The goal of the experiment is to compare the mean breaking strength of fibers produced by the three machines. When testing this, the technician also measured the diameter of each fiber, and this measurement may be used as a concomitant variable to improve precision of the estimates.
}
\source{
Montgomery, D. C. (2013) \emph{Design and Analysis of Experiments} (8th ed.). John Wiley and Sons, ISBN 978-1-118-14692-7.
}

\examples{
require(lsmeans)
fiber.lm <- lm(strength ~ diameter + machine, data=fiber)
ref.grid(fiber.lm)

# Covariate-adjusted means and comparisons
lsmeans(fiber.lm, pairwise ~ machine)
}
\keyword{datasets}
