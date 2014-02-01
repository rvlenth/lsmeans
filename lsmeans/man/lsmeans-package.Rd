\name{lsmeans-package}
\alias{lsmeans-package}
\docType{package}
\title{
Least-squares means
}
\description{
This package provides methods for obtaining so-called least-squares means for factor combinations in a variety of fitted linear models. It can also compute contrasts or linear combinations of these least-squares means, (several standard contrast families are provided), and in addition can estimate and contrasts slopes of trend lines.
Popular adjustments for multiple-comparisons are provided. There is also an interface to \code{\link[multcomp]{glht}}, if the \pkg{multcomp} package is installed, and a compact letter display method for pairwise comparisons.

Examples and discussion are available via \code{vignette("using-lsmeans", package="lsmeans")}.

Starting with \pkg{lsmeans} version 2, a new object framework based on \emph{reference grids} is used that increases flexibility and provides for extending its capabilities to additional model objects. Use \code{vignette("lsmeans-changes")} for information on the user impact of these changes.
See \link{extending-lsmeans} for details on extending to other models.
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
\link[=ref.grid-class]{The ref.grid class},
\code{\link[multcomp:glht]{multcomp}}
}
