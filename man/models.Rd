%    Copyright (c) 2012-2016 Russell V. Lenth                                %

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
  \item{mlm, maov, manova}{When there is a multivariate response, the different responses are treated as if they were levels of a factor -- named \code{rep.meas} by default. The \code{mult.name} argument may be used to change this name. The \code{mult.levs} argument may specify a named list of one or more sets of levels. If this has more than one element, then the multivariate levels are expressed as combinations of the named factor levels via the function \code{\link{expand.grid}}.}
  \item{aovlist}{Support for these objects is limited. To avoid strong biases in the predictions, the \code{contrasts} attribute of all factors should be of a type that sums to zero -- for example, \code{"contr.sum"}, \code{"contr.poly"}, or \code{"contr.helmert"} but \emph{not} \code{"contr.treatment"}.  Only intra-block estimates of covariances are used. That is, if a factor appears in more than one error stratum, only the covariance structure from its lowest stratum is used in estimating standard errors. Degrees of freedom are obtained using the Satterthwaite method. In general, \code{aovList} support is best with balanced designs, and due caution in the use of contrasts. If a \code{vcov.} argument is supplied, it must yield a single covariance matrix for the unique fixed effects, and the degrees of freedom are set to \code{NA}.}
}} %stats


%%% !!! ALPHABETICAL FROM HERE !!!

\section{\pkg{afex} package}{
\describe{
\item{mixed}{Support for \code{mixed} objects has been removed. Version 0.14 and later of \pkg{afex} provides new object classes with their own \pkg{lsmeans} support.}}}
% afex
%%%Support for the \code{full.model} element of these objects is the same as that for \code{merMod} in the \pkg{lme4} package -- see below. However, for \pkg{afex} versions 0.10-113 and earlier, the \code{data} argument is required in calls to \code{lsmeans} or \code{ref.grid}, as the information about the original dataset is not preserved in the object.}


\section{\pkg{betareg} package}{
\describe{
\item{betareg}{The additional \code{mode} argument has possible values of \code{"response"}, \code{"link"}, \code{"precision"}, \code{"phi.link"}, \code{"variance"}, and \code{"quantile"}, which have the same meaning as the \code{type} argument in \code{predict.betareg} -- with the addition that \code{"phi.link"} is like \code{"link"}, but for the precision portion of the model. When \code{mode = "quantile"} is specified, the additional argument \code{quantile} (a numeric scalar or vector) specifies which quantile(s) to compute; the default is 0.5 (the median). Also in \code{"quantile"} mode, an additional variable \code{quantile} is added to the reference grid, and its levels are the values supplied.
}
}}% betareg


\section{\pkg{CARBayes} package}{
\describe{
\item{carbayes}{The user \emph{must} supply (via the \code{data} argument) the dataset used in fitting the model. As with other MCMC-based objects, the summaries and such are frequentist, but the \code{as.mcmc} method provides a posterior sample of the desired quantities.}
}}% CARBayes


\section{\pkg{coxme} package}{
\describe{
\item{coxme}{No additional options. Support for these models is experimental; may throw errors or incorrect results.}
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
}
If a matrix or function is supplied as \code{vcov.method}, it is interpreted as a \code{vcov.} specification as described for \code{...} in \code{\link{ref.grid}}.
}% geepack


\section{\pkg{glmmADMB} package}{
\describe{
\item{glmmadmb}{No extended features.}
}}% glmmadmb



\section{\pkg{lme4} package}{
\describe{
\item{lmerMod}{There is an optional \code{mode} argument that defaults to \code{get.lsm.option("lmer.df")} (which in turn defaults to \code{"satterthwaite"}). The possible values are "satterthwaite", "kenward-roger", and "asymptotic" (these are partially matched and case-insensitive). With \code{"satterthwaite"}, d.f. are obtained using code from the \pkg{lmerTest} package, if installed. With \code{"kenward-roger"}, d.f. are obtained using code from the \pkg{pbkrtest} package, if installed. With \code{"asymptotic"}, or if the needed package is not installed, d.f. are set to \code{NA}.

A by-product of the Kenward-Roger method is that the covariance matrix is adjusted using \code{\link[pbkrtest]{vcovAdj}}. This can require considerable computation; so to avoid that overhead, the user should opt for the Satterthwaite or asymptotic method; or, for backward compatibility, may disable the use of \pkg{pbkrtest} via \samp{lsm.options(disable.pbkrtest=TRUE)} (this does not disable the \pkg{pbkrtest} package entirely, just its use in \pkg{lsmeans}). The computation time required depends roughly on the number of observations, \emph{N}, in the design matrix (because a major part of the computation involves inverting an \emph{N x N} matrix). Thus, \pkg{pbkrtest} is automatically disabled if \emph{N} exceeds the value of \code{get.lsm.option("pbkrtest.limit")}. If desired, the user may use \code{lsm.options} to adjust this limit from the default of 3000.

The \code{df} argument may be used to specify some other degrees of freedom. Note that if \code{df} and \code{method = "satterthwaite"} are both specified, the covariance matrix is adjusted but the K-R degrees of freedom are not used.
}
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

%\section{\pkg{mgcv} package}{
%\describe{
%\item{gam}{Supported by virtue of inheritance from \code{glm}.}
%\item{gamm}{Not supported at this time.}
%}}% mgcv

\section{\pkg{MCMCglmm} package}{
\describe{
\item{MCMCglmm}{Currently, I have found no way to reconstruct the data based on information in the object; thus, you \emph{must} provide the dataset via the \code{data} argument. In addition, the \code{contrasts} specifications are not recoverable from the object, so the system default must match what was actually used in fitting the model. The usual \code{summary}, \code{test}, etc. methods provide frequentist analyses of the results based on the posterior means and covariances. However, an \code{as.mcmc} method is provided that creates an \code{mcmc} object that can be summarized or plotted using the \pkg{coda} package. It provides a posterior sample of lsmeans for the given reference grid, based on the posterior sample of the fixed effects from the \code{MCMCglmm} object.}
}}% MCMCglmm

\section{\pkg{MCMCpack} package (and perhaps others)}{
\describe{
\item{mcmc}{Certain linear-model or mixed-model objects are of class \code{mcmc}, and contain a sample from the posterior distribution of fixed-effect coefficients. In some cases (e.g., results of \code{MCMCregress} and \code{MCMCpoisson}), the object may include a \code{"call"} attribute that \code{lsmeans} can use to reconstruct the data and obtain a basis for the least-squares means. If not, a \code{formula} and \code{data} argument are provided that may help produce the right results. In addition, the \code{contrasts} specifications are not recoverable from the object, so the system default must match what was actually used in fitting the model. As for other MCMC-based objects, the summaries and such are frequentist, but the \code{as.mcmc} method provides a posterior sample of the desired quantities.}
}}% MCMCpack


\section{\pkg{nlme} package}{
\describe{
\item{gls}{No additional features. Degrees of freedom are computed using \code{N - p} in \code{object$dims}. This is consistent with \code{nlme:::summary.gls} but seems questionable.}
\item{lme}{Degrees of freedom are obtained using a containment method, i.e., the minimum of those elements of \code{object$fixDF$X} receiving nonzero weight (but with a correction to the \code{lme} object's intercept df). (This is similar to \pkg{SAS}'s containment method, but I believe \pkg{SAS} does it incorrectly when the estimands are not contrasts.) The optional argument \code{sigmaAdjust} (defaults to \code{TRUE}) will adjust standard errors like in \code{\link[nlme]{summary.lme}} when the model is fitted using the \code{"ML"} method. \bold{Note:} \code{sigmaAdjust} is comparable to \code{adjustSigma} in \code{\link[nlme]{summary.lme}} but it is renamed to avoid conflicting with \code{adjust}.}
\item{nlme}{Support is provided for inferences on parameters named in the \code{fixed} part of the model. The user \emph{must} specify \code{param} in the call and give the name of a parameter that appears in the right-hand side of a \code{fixed} formula. Degrees of freedom are obtained using the containment-like method described above for \code{lme}.}
}}% nlme

\section{\pkg{nnet} package}{
\describe{
\item{multinom}{
The reference grid includes a pseudo-factor with the same name and levels as the multinomial response. There is an optional \code{mode} argument which should match \code{"prob"} or \code{"latent"}. With \code{mode = "prob"}, the reference-grid predictions consist of the estimated multinomial probabilities. The \code{"latent"} mode returns the linear predictor, recentered so that it averages to zero over the levels of the response variable (similar to sum-to-zero contrasts). Thus each latent variable can be regarded as the log probability at that level minus the average log probability over all levels. 

Please note that, because the probabilities sum to 1 (and the latent values sum to 0) over the multivariate-response levels, all sensible results from \code{lsmeans} must involve that response as one of the factors. For example, if \code{resp} is a response with \eqn{k} levels, \code{lsmeans(model, ~ resp | trt)} will yield the estimated multinomial distribution for each \code{trt}; but \code{lsmeans(model, ~ trt)} will just yield the average probability of \eqn{1/k} for each \code{trt}.
}}}% nnet, multinom


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


\section{\pkg{pscl} package}{
\describe{
\item{hurdle, zeroinfl}{
Two optional arguments -- \code{mode} and \code{lin.pred} -- are provided. The \code{mode} argument has possible values \code{"response"} (the default), \code{"count"}, \code{"zero"}, or \code{"prob0"}. \code{lin.pred} is logical and defaults to \code{FALSE}.

With \code{lin.pred = FALSE}, the results are comparable to those returned by \code{predict(..., type = "response")}, \code{predict(..., type = "count")}, \code{predict(..., type = "zero")}, or \code{predict(..., type = "prob")[, 1]}. See the documentation for \code{\link[pscl]{predict.hurdle}} and \code{\link[pscl]{predict.zeroinfl}}.

The option \code{lin.pred = TRUE} only applies to \code{mode = "count"} and \code{mode = "zero"}. The results returned are on the linear-predictor scale, with the same transformation as the link function in that part of the model. The predictions for a reference grid with \code{mode = "count"}, \code{lin.pred = TRUE}, and \code{type = "response"} will be the same as those obtained with \code{lin.pred = FALSE} and \code{mode = "count"}; however, any LS means derived from these grids will be defferent, because the averaging is done on the log-count scale and the actual count scale, respectively -- thereby producing geometric means versus arithmetic means of the predictions.

If the \code{vcov.} argument is used (see details in \code{\link{ref.grid}}), it must yield a matrix of the same size as would be obtained using \code{\link[pscl]{vcov.hurdle}} or \code{\link[pscl]{vcov.zeroinfl}} with its \code{model} argument set to \code{("full", "count", "zero")} in respective correspondence with \code{mode} of \code{("mean", "count", "zero")}.
If \code{vcov.} is a function, it must support the \code{model} argument.
}
}}% pscl


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
