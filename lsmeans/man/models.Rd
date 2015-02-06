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


\section{\pkg{gam} package}{
\describe{
\item{gam}{Currently, \code{gam} objects are not supported. Past versions of \pkg{lsmeans} appeared to support \code{gam} models owing to inheritance from \code{lm}, but the results were incorrect because spline features were ignored. We now explicitly trap \code{gam} objects to avoid tyhese misleading analyses.}
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

\section{\pkg{lme4.0} package}{
\describe{
\item{mer}{Only asymptotic results are available (no d.f.).}
}}% lme4.0

\section{\pkg{MASS} package}{
\describe{
\item{glmmPQL}{Supported by virtue of inheritance from \code{lme} in the \pkg{nlme} package.}
\item{glm.nb}{Supported by virtue of inheritance from \code{glm}.}
\item{polr}{There are two optional arguments: \code{mode} and \code{rescale} (which defaults to \samp{c(0,1)}). For details, see the documentation below regarding the support for the \pkg{ordinal} package, which produces comparable objects (but since \code{polr} does not support scale models, \code{mode="scale"} is not supported).
Tests and confidence intervals are asymptotic.}
\item{rlm}{Supported by virtue of inheritance from \code{lm}.}
}}% MASS

\section{\pkg{mgcv} package}{
\describe{
\item{gam}{Supported by virtue of inheritance from \code{glm}.}
\item{gamm}{Not supported at this time.}
}}% mgcv


\section{\pkg{nlme} package}{
\describe{
\item{gls}{No degrees of freedom are available for these objects, so tests and confidence intervals are asymptotic.}
\item{lme}{No degrees of freedom, as per \code{gls} objects. The optional argument \code{adjustSigma} (defaults to \code{TRUE}) will adjust standard errors like in \code{\link[nlme]{summary.lme}} when the model is fitted using the \code{"ML"} method. \bold{Note:} It is possible for the \code{adjust} argument (for \eqn{p}-value adjustments) to conflict with \code{adjustSigma}. The workaround is to specify both: e.g., \samp{lsmeans(mod.lme, pairwise ~ trt, adjustSigma = TRUE, adjust = "none")}}
}}% nlme

\section{\pkg{ordinal} package}{
\describe{
\item{clm,clmm}{The reference grid will include all variables that appear in the main model as well as those in the \code{scale} or \code{nominal} models. There are two optional arguments: \code{mode} (a character string) and \code{rescale} (which defaults to \samp{c(0,1)}). \code{mode} should match one of \code{"latent"} (the default), \code{"linear.predictor"}, \code{"cum.prob"}, \code{"exc.prob"}, \code{"prob"}, \code{"mean.class"}, or \code{"scale"}. 

With \samp{mode = "latent"}, the reference-grid predictions are made on the scale of the latent variable implied by the model. The scale and location of this latent variable are arbitrary, and may be altered via \code{rescale}. The predictions are multiplied by \samp{rescale[2]}, then \samp{rescale[1]} is added. Keep in mind that the scaling is related to the link function used in the model; for example, changing from a probit link to a logistic link will inflate the latent values by around \eqn{\pi/\sqrt{3}}{pi/sqrt(3)}, all other things being equal. \code{rescale} has no effect for other values of \code{mode}.

With \samp{mode = "linear.predictor"} \code{mode = "cum.prob"}, and \code{mode = "exc.prob"}, the boundaries between categories (i.e., thresholds) in the ordinal response are included in  the reference grid as a pseudo-factor named \code{cut}. The reference-grid predictions are then of the cumulative probabilities at each threshold (for \code{mode = "cum.prob"}), exceedance probabilities (one minus cumulative probabilities, for \code{mode = "exc.prob"}), or the link function thereof (for \code{mode = "linear.predictor"}).

With \code{mode = "prob"}, a pseudo-factor with the same name as the model's response variable is created, and the grid predictions are of the probabilities of each class of the ordinal response. With \code{"mean.class"}, the returned results are means of the ordinal response, interpreted as a numeric value from 1 to the number of classes, using the \code{"prob"} results as the estimated probability distribution for each case.

With \code{mode = "scale"}, and the fitted object incorporates a scale model, least-squares means are obtained for the factors in the scale model instead of the response model. The grid is constructed using only the factors in the scale model.

Any grid point that is non-estimable by either the location or the scale model (if present) is set to \code{NA}, and any LS-means involving such a grid point will also be non-estimable. A consequence of this is that if there is a rank-deficient \code{scale} model, and then \emph{all} latent responses become non-estimable because the predictions are made using the average log-scale estimate.

Tests and confidence intervals are asymptotic.}
}}% ordinal

\section{\pkg{rms} package}{
\describe{
\item{Potential masking issue}{
Both \pkg{rms} and \pkg{lsmeans} offer \code{contrast} methods, and whichever package is loaded later masks the other. Thus, you may need to call \code{lsmeans::contrast} or \code{rms::contrast} explicitly to access the one you want. 
}
\item{Objects inheriting from rms}{
Standard support is provided. However, with models having more than one intercept (e.g. from \code{orm}), a \code{mode} argument is provided that works similarly to that for the \pkg{ordinal} package. The available modes are \code{"middle"} (the default), \code{"latent"}, \code{"linear.predictor"}, \code{"cum.prob"}, \code{"exc.prob"}, \code{"prob"}, and \code{"mean.class"}. All are as described for the \pkg{ordinal} package, except as noted below.

With \code{mode = "middle"} (this is the default), the middle intercept is used, comparable to the default for \pkg{rms}'s \code{Predict} function. This is quite similar in concept to \code{mode = "latent"}, where all intercepts are averaged together.

Results for \code{mode = "linear.predictor"} are reversed from those in the \pkg{ordinal} package, because \code{orm} models predict the link function of the \emph{upper}-tail (exceedance) probabilities.

With \code{mode = "prob"}, a pseudo-factor is created having the same name as the model response variable, but its levels are always integers \samp{1, 2, ...} regardless of the levels of the original response.
}
}}% rms


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
