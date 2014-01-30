\name{lsmeans-package}
\alias{lsmeans-package}
\docType{package}
\title{
Least-squares means
}
\description{
This package provides methods for obtaining so-called least-squares means for factor combinations in a variety of fitted linear models. It can also compute contrasts or linear combinations of these least-squares means, (several standard contrast families are provided), and in addition can estimate and contrasts slopes of trend lines.
Popular adjustments for multiple-comparisons are provided. There is also an interface to \code{\link[multcomp]{glht}}, if the \pkg{multcomp} package is installed, and a coompact letter display method for pairwise comparisons.

In analyzing models fitted by the \pkg{lme4} package, if the \code{pbkrtest} package is installed, its \code{\link[pbkrtest]{vcovAdj}} function is used to bias-adjust the covariance matrix and obtain degrees of freedom.

Examples and discussion are available via \code{vignette("using-lsmeans", package="lsmeans")}.

It is possible to extend the package's capabilities to additional types of models by writing methods that provide the needed information. See \link{extending-lsmeans} for details.
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
\code{\link[lme4:lmer]{lme4}}, 
\code{\link[multcomp:glht]{multcomp}}, 
\code{\link[pbkrtest:vcovAdj]{pbkrtest}}
}
