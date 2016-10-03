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
}

\section{Overview}{
\describe{
\item{Concept}{
Least-squares means (see Searle \emph{et al.} 1980, who prefer the term \dQuote{predicted marginal means} (PMM)) are popular for summarizing linear models that include factors. For balanced experimental designs, they are just the marginal means. For unbalanced data, they in essence estimate what you \emph{would} have observed that the data arisen from a  balanced experiment.
}
\item{Reference grids}{
The implementation in \pkg{lsmeans} relies on our own concept of a \emph{reference grid}, which is an array of factor and predictor levels.  
Predictions are made on this grid, and least-squares means are defined as averages of these predictions over zero or more dimensions of the grid. The function \code{\link{ref.grid}} explicitly creates a reference grid (\code{ref.grid} object) that can subsequently be used to obtain least-squares means. The \code{\link{update}} method is used to change its properties.

Our reference-grid framework expands slightly upon Searle \emph{et al.}'s definitions of PMMs, in that it is possible to include multiple levels of covariates in the grid.
}

\item{Models supported}{
Many linear models are supported by the package, including \code{lm}, \code{glm}, \code{aovList}, and \code{mlm} in the \pkg{stats} package, as well as fitted-model objects from several contributed packages including \pkg{nlme}, \pkg{lme4}, \pkg{survival}, and \pkg{geepack}. The help page for \code{\link{models}} provides more details, including, in some cases, additional \code{ref.grid} arguments that might affect the subsequent analysis. Also, some models require other packages be installed in order to  obtain all the available features.
}
\item{Least-squares means}{
The \code{\link{lsmeans}} function computes least-squares means given a \code{ref.grid} object or a fitted model, and a specification indicating what factors to include. The \code{\link{lstrends}} function creates the same sort of results for estimating and comparing slopes of fitted lines. Both return an \code{lsmobj} object very much like a reference grid, but with possibly fewer factors involved.
}
\item{Summaries and analysis}{
The \code{\link{summary}}  method may be used to display a \code{ref.grid} or an \code{lsmobj}. Special-purpose summaries are available via \code{\link{confint}} and \code{\link{test}}, the latter of which can also do a joint test of several estimates. The user may specify by variables, multiplicity-adjustment methods, confidence levels, etc., and if a transformation or link function is involved, may reverse-transform the results to the response scale.
}
\item{Contrasts and comparisons}{
The \code{\link{contrast}} method is used to obtain contrasts among the estimates; several standard contrast families are available such as deviations from the mean, polynomial contrasts, and comparisons with one or more controls. Another \code{lsmobj} object is returned, which can be summarized or further analyzed. For convenience, a \code{pairs} method is provided for the case of pairwise comparisons. Related to this is the \code{\link{cld}} method, which provides a compact letter display for grouping pairs of means that are not significantly different. \code{cld} requires the \pkg{multcompView} package.
}
\item{Graphs}{The \code{\link[=plot.lsmobj]{plot}} method will display side-by-side confidence intervals for the estimates, and/or \sQuote{comparison arrows} whereby the significance of pairwise differences can be judged by how much they overlap. The \code{\link{lsmip}} function displays estimates like an interaction plot, multi-paneled if there are by variables. These graphics capabilities require the \pkg{lattice} package be installed.
}
\item{\pkg{multcomp} interface}{
The \code{\link{as.glht}} function and \code{\link{glht}} method for \code{lsmobj}s provide an interface to the \code{\link[multcomp]{glht}} function in the \pkg{multcomp} package, thus providing for more exacting simultaneous estimation or testing. The package also provides an \code{\link{lsm}} method that works as an alternative to \code{\link[multcomp]{mcp}} in a call to \code{glht}.
}
}}% overview

\section{Additional information}{
Examples and discussion are available via \code{vignette("using-lsmeans", package="lsmeans")}.

Some features of the \pkg{lsmeans} require (or are enhanced by) additional packages that are loaded when needed. Since they are not \dQuote{required} packages, they are not automatically installed with \pkg{lsmeans}. We highly recommend that users also install the following packages: \pkg{multcomp} (if \code{\link{cld}}, \code{\link{glht}}, or \code{\link{as.glht}} are to be used), \pkg{multcompView} (for \code{\link{cld}}), \pkg{lattice} (for \code{\link{plot}} and \code{\link{lsmip}}), and \pkg{lmerTest} or \pkg{pbkrtest} (for models fitted by the \pkg{lme4} package).

Starting with \pkg{lsmeans} version 2, a new object framework based on \emph{reference grids} is used that increases flexibility and provides for extending its capabilities to additional model objects. Use \code{vignette("lsmeans-changes")} for information on the user impact of these changes.

It is possible to write your own interfaces for models not yet supported by \pkg{lsmeans}. See the help page \code{\link{extending-lsmeans}} and \code{vignette("extending")} for details on how to do this.
} % add'l information



\details{
\tabular{ll}{
Package: \tab lsmeans\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
Other information: \tab See DESCRIPTION\cr
}
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
