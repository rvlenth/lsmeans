\name{recover.data}
\alias{recover.data}
\alias{recover.data.call}
\alias{lsm.basis}
%%\alias{.all.vars} % documented in the vignette
%%\alias{nonest.basis} % moved to estimability package
%%\alias{is.estble}

% additional ref topic
\alias{extending-lsmeans}


\title{Support functions for creating a reference grid}
\description{
This documents the methods used to create a \code{\link[=ref.grid-class]{ref.grid}} object from a fitted model.
}
\usage{
recover.data(object, ...)
\S3method{recover.data}{call}(object, trms, na.action, 
    data = NULL, params = NULL, ...)

lsm.basis(object, trms, xlev, grid, ...)
%%%
%%%.all.vars(expr, retain = c("\\\\$", "\\\\[\\\\[", "\\\\]\\\\]"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object returned from a model-fitting function.}
  \item{trms}{The \code{\link{terms}} component of \code{object}}
  \item{xlev}{Named list of levels of factors in the model frame. This should \emph{not} include levels of factors created in the model itself, e.g., by including a \code{factor} call in the model formula.}
  \item{grid}{A \code{data.frame} containing predictor values at which predictions are needed.}
  \item{na.action}{Integer vector of indices of observations to ignore; or \code{NULL} if none}
  \item{data}{Data frame. Usually, this is \code{NULL}. However, if non-null, this is used in place of the reconstructed dataset. It must have all of the predictors used in the model, and any factor levels must match those used in fitting the model.}
  \item{params}{Character vector giving the names of any variables in the
  model formula that are \emph{not} predictors. An example would be a   
  variable \code{knots} specifying the knots to use in a spline model.
}
%  \item{expr}{A formula}
%  \item{retain}{Character vector of operators to retain (escaped as for \code{\link{gsub}})}
  \item{\dots}{Additional arguments passed to other methods.}
}

\details{
To create a reference grid, the \code{ref.grid} function needs to reconstruct the data used in fitting the model, and then obtain a matrix of linear functions of the regression coefficients for a given grid of predictor values. These tasks are performed by calls to \code{recover.data} and \code{lsm.basis} respectively. 

To extend \pkg{lsmeans}'s support to additional model types, one need only write S3 methods for these two functions. The existing methods serve as helpful guidance for writing new ones.  Most of the work for \code{recover.data} can be done by its method for class \code{"call"}, providing the \code{terms} component and \code{na.action} data as additional arguments. Writing an \code{lsm.basis} method is more involved, but the existing methods (e.g., \code{lsmeans:::lsm.basis.lm}) can serve as models. See the ``Value'' section below for details on what it needs to return. Also, certain \code{recover.data} and \code{lsm.basis} methods are exported from \pkg{lsmeans}, so if your object is based on another model-fitting object, it may be that all that is needed is to call one of these exported methods and perhaps make modifications to the results. Contact the developer if you need others of these exported.

If the model has a multivariate response, \code{bhat} needs to be \dQuote{flattened} into a single vector, and \code{X} and \code{V} must be constructed consistently. 

In models where a non-full-rank result is possible (often you can tell by seeing if there is a \code{singular.ok} argument in the model-fitting function), \code{summary} and \code{predict} check the estimability of each prediction, using the \code{\link[estimability]{nonest.basis}} function in the \pkg{estimability} package. 

The models already supported are detailed in \code{\link{models}}. Some packages may provide additional \pkg{lsmeans} support for its object classes. 
}


\value{
\code{recover.data} should return a \code{data.frame} containing all the variables in the original data that appear as predictors in the model. Several attributes need to be included as well; see the code for \code{lsmeans:::recover.data.lm}.

\code{lsm.basis} should return a \code{list} with the following elements:
\item{X}{The matrix of linear fcns over \code{grid}, having the same number of rows as \code{grid} and the number of columns equal to the length of \code{bhat}.}
\item{bhat}{The vector of regression coefficients for fixed effects. This should \emph{include} any \code{NA}s that result from rank deficiencies.}
\item{nbasis}{A matrix whose columns form a basis for non-estimable functions of beta, or a 1x1 matrix of \code{NA} if there is no rank deficiency.}
\item{V}{The estimated covariance matrix of \code{bhat}.}
\item{dffun}{A function of \code{(k, dfargs)} that returns the degrees of freedom associated with \code{sum(k * bhat)}.}
\item{dfargs}{A \code{list} containing additional arguments needed for \code{dffun}.}

%%%\code{.all.vars} is an enhancement of \code{\link{all.vars}}, whereby the operators specified in \code{retain} are left intact. Thus, \code{All.vars(foo$y ~ bar[[2]])} returns \code{"foo$y", "bar[[2]]"}, whereas \code{all.vars} returns \code{"foo", "y", "bar"}
}

%\references{}

\section{Optional hooks}{
Some models may need something other than standard linear estimates and standard errors. If so, custom functions may be pointed to via the items \code{misc$estHook}, \code{misc$vcovHook} and \code{misc$postGridHook}. If just the name of the hook function is provided as a character string, then it is retrieved using \code{\link{get}}.

The \code{estHook} function should have arguments \samp{(object, do.se, tol, ...)} where \code{object} is the \code{ref.grid} or \code{lsmobj} object, \code{do.se} is a logical flag for whether to return the standard error, and \code{tol} is the tolerance for assessing estimability. It should return a matrix with 3 columns: the estimates, standard errors (\code{NA} when \code{do.se==FALSE}), and degrees of freedom (\code{NA} for asymptotic). The number of rows should be the same as \samp{object@linfct}. The \code{vcovHook} function should have arguments \samp{(object, tol, ...)} as described. It should return the covariance matrix for the estimates. Finally, \code{postGridHook}, if present, is called at the very end of \code{ref.grid}; it takes one argument, the constructed \code{object}, and should return a suitably modified\code{ref.grid} object.
}

\section{Additional functions}{
A few additional functions used in the \pkg{lsmeans} codebase are exported as they may be useful to package developers. See details near the end of the vignette \code{"extending"}.
}

\author{
Russell V. Lenth
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{models}}, \code{\link{ref.grid}}, \code{\link{ref.grid-class}}
}
\examples{
\dontrun{
    require(lsmeans)

    # Fit a 2-factor model with two empty cells
    warpsing.lm <- lm(breaks ~ wool*tension, 
        data = warpbreaks, subset = -(16:40))

    lsmeans:::recover.data.lm(warpsing.lm, data = NULL)
    grid = with(warpbreaks, 
        expand.grid(wool = levels(wool), tension = levels(tension)))
    lsmeans:::lsm.basis.lm(warpsing.lm, delete.response(terms(warpsing.lm)),
        warpsing.lm$xlevels, grid)
} % end dontrun

}


\keyword{ models }
\keyword{ regression }
