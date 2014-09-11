\name{models}
\alias{models}

\title{Models supported in \pkg{lsmeans}}

\description{
Here we document what model objects may be used with \pkg{lsmeans}, and some special features of some of them. We start with those in the \pkg{stats} package; the other packages follow in alphabetical order.

Certain objects are affected by optional arguments to functions that construct \code{ref.grid} or \code{lsmobj} objects, including \code{\link{ref.grid}}, \code{\link{lsmeans}}, \code{\link{lstrends}}, and \code{\link{lsmip}}. When \dQuote{arguments} are mentioned in the subsequent object-by-object documentation, we are talking about arguments in these constructors.

Additional models can be supported by writing appropriate \code{recover.data} and \code{lsm.basis} methods. See \code{\link{extending-lsmeans}} and \code{vignette("extending")} for details.
}


\section{\pkg{stats} package}{
  \describe{
  \item{lm, aov, glm}{No extended features. Note that the \code{lm} support often extends to a number of model objects that inherit from it, such as \code{rlm} in the \pkg{MASS} package and \code{rsm} in the \pkg{rsm} package.}
  \item{mlm}{When there is a multivariate response, the different responses are treated as if they were levels of a factor -- named \code{rep.meas} by default. The \code{mult.name} argument may be used to change this name. The \code{mult.levs} argument may specify a named list of one or more sets of levels. If this has more than one element, then the multivariate levels are expressed as combinations of the named factor levels via the function \code{\link{expand.grid}}.}
  \item{aovList}{Support for these objects is limited. To avoid strong biases in the predictions, the \code{contrasts} attribute of all factors should be of a type that sums to zero -- for example, \code{"contr.sum"}, \code{"contr.poly"}, or \code{"contr.helmert"} but \emph{not} \code{"contr.treatment"}.  Only intra-block estimates of covariances are used. That is, if a factor appears in more than one error stratum, only the covariance structure from its lowest stratum is used in estimating standard errors. In general, \code{aovList} support is best with balanced designs, and due caution in the use of contrasts.}
}} %stats


%%% !!! ALPHABETICAL FROM HERE !!!

\section{\pkg{afex} package}{
\describe{
\item{mixed}{Support for the \code{full.model} element of these objects is the same as that for \code{merMod} in the \pkg{lme4} package -- see below. However, for \pkg{afex} versions 0.10-113 and earlier, the \code{data} argument is required in calls to \code{lsmeans} or \code{ref.grid}, as the information about the original dataset is not preserved in the object.}
}}% afex


\section{\pkg{coxme} package}{
\describe{
\item{coxme}{No extended features.}
}}% coxme



\section{\pkg{gee} and \pkg{geepack} packages}{
These models all have more than one covariance estimate available, and it may be selected by supplying a string as the \code{vcov.method} argument. It is partially matched with the available choices; thus, for example, \samp{vcov = "n"} translates to \samp{vcov.method = "naive"}
\describe{
\item{gee}{Available covariance estimates are specified in \code{vcov.method} as \code{"robust"} (the default) and \code{"naive"}.}
\item{geeglm, geese}{Available covariance estimates are specified in \code{vcov.method} as \code{"vbeta"} (the default), \code{"vbeta.naiv"}, \code{"vbeta.j1s"}, or \code{"vbeta.fij"}. The aliases \code{"robust"} (for \code{"vbeta"}) and \code{"naive"} (for \code{"vbeta.naiv"} are also accepted.}
}}% geepack


\section{\pkg{glmmADMB} package}{
\describe{
\item{glmmadmb}{No extended features.}
}}% glmmadmb



\section{\pkg{lme4} package}{
\describe{
\item{lmerMod}{If the \pkg{pbkrtest} package is installed, degrees of freedom for confidence intervals and tests are obtained using its \code{\link[pbkrtest]{ddf_Lb}} function, and the covariance matrix is adjusted using \code{\link[pbkrtest]{vcovAdj}}.
If \pkg{pbkrtest} is not installed, the covariance matrix is not adjusted, degrees of freedom are set to \code{NA}, and asymptotic results are displayed. 

The user may disable the use of \pkg{pbkrtest} via \samp{lsm.options(disable.pbkrtest=TRUE)} (this does not disable the \pkg{pbkrtest} package entirely, just its use in \pkg{lsmeans}). The \code{df} argument may be used to specify some other degrees of freedom. Specifying \code{df} is \emph{not} equivalent to disabling \pkg{pbkrtest}, because if not disabled, the covariance matrix is still adjusted.}
\item{glmerMod}{No degrees of freedom are available for these objects, so tests and confidence intervals are asymptotic.}
}}% lme4


\section{\pkg{MASS} package}{
\describe{
\item{polr}{The boundaries between categories in the ordinal response are included in  the reference grid as the factor \code{cut}. Thus, the user may specify which cut-point to use in estimating cumulative probabilities.}
\item{rlm}{Supported by virtue of inheritance from \code{lm}.}
}}% MASS



\section{\pkg{nlme} package}{
\describe{
\item{lme, gls}{No degrees of freedom are available for these objects, so tests and confidence intervals are asymptotic.}
}}% nlme


\section{\pkg{survival} package}{
\describe{
\item{survreg, coxph}{No extended features.}
}}% survival



%\references{}
\author{
Russell V. Lenth
}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ref.grid}}, \code{\link{lsm.basis}}
}

\keyword{ models }
\keyword{ regression }
\keyword{ htest }
