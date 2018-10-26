%    Copyright (c) 2012-2016 Russell V. Lenth                                %

\name{lsmeans-package}
\alias{lsmeans-package}
\docType{package}
\encoding{utf-8}
\title{
Least-squares means
}
\description{
This package provides methods for obtaining so-called least-squares means for factor combinations in a variety of fitted linear models. It can also compute contrasts or linear combinations of these least-squares means, (several standard contrast families are provided), and in addition can estimate and contrast slopes of trend lines.
Popular adjustments for multiple-comparisons are provided, as well as graphical ways of displaying the results.

Almost the entire codebase for \pkg{lsmeans} now resides in the \pkg{emmeans} package 
(named for the more general term, \dQuote{estimated marginal means}). \pkg{lsmeans} exists only as a transitional entity for the few remaining packages that depend on it. 
}

\references{
Russell V. Lenth (2016)
Least-Squares Means: The R Package lsmeans.
\emph{Journal of Statistical Software}, 69(1), 1-33. doi:10.18637/jss.v069.i01
  
Searle S.R. Speed F.M. Milliken G.A. (1980) 
Population marginal means in the linear model: An alternative to least squares means. 
\emph{The American Statistician} \bold{34}(4), 216-221.
}

\author{
Russell V. Lenth (author), Maxime Hervé (contributor)

Maintainer: Russ Lenth <russell-lenth@uiowa.edu>
}


\keyword{ package }
\keyword{ models }
\keyword{ regression }
\keyword{ htest }
