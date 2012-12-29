\name{lsmeans-package}
\alias{lsmeans-package}
\docType{package}
\title{
Least-squares means
}
\description{
This package provides the function \code{\link{lsmeans}} for outputting so-called least-squares means for factor combinations in an \code{lm}, \code{lme}, \code{gls}, or \code{lmer} object. It can also compute contrasts or linear combinations of these least-squares means, and several standard contrast families are provided.
It may also be used for the obtaining linear predictors and contrasts thereof in \code{glm} and \code{glmer} objects.
Popular adjustments for multiple-comparisons are provided; optionally, contrasts may instead be passed to \code{\link[multcomp]{glht}}, if the \code{multcomp} package is installed. 

In analyzing \code{mer} objects, if the \code{pbkrtest} package is installed, its \code{\link[pbkrtest]{vcovAdj}} function is used to bias-adjust the covariance matrix and obtain degrees of freedom.

Examples and discussion are available via \code{vignette("using-lsmeans", package="lsmeans")}.
}
\details{
\tabular{ll}{
Package: \tab lsmeans\cr
Type: \tab Package\cr
Version: \tab See DESCRIPTION\cr
Date: \tab See DESCRIPTION\cr
License: \tab GPL-2\cr
}
See the documentation for \code{\link{lsmeans}}.
}
\author{
Russell V. Lenth

Maintainer: Russ Lenth <russell-lenth@uiowa.edu>
}


\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ htest }
\seealso{
\code{\link[nlme:lme]{nlme}}
\code{\link[lme4:lmer]{lme4}}
\code{\link[multcomp:glht]{multcomp}}
\code{\link[pbkrtest:vcovAdj]{pbkrtest}}
}
