%    Copyright (c) 2012-2016 Russell V. Lenth                                %

\name{lsmeans}
\alias{lsmeans}
\alias{lsmeans.formula}
\alias{lsmeans.list}
\alias{lsmeans.character}
\alias{lsmeans.character.ref.grid}
\alias{lstrends}
\alias{lsmobj}
\alias{pmmeans}
\alias{pmtrends}
\alias{pmmobj}

\title{Least-squares means (or predicted marginal means)}
\description{
Compute least-squares means (predicted marginal means) for specified factors or factor combinations in a linear model,
and optionally comparisons or contrasts among them.
}
\usage{
\method{lsmeans}{character}(object, specs, ...)
## (used when 'specs' is 'character')

\method{lsmeans}{character.ref.grid}(object, specs, by = NULL, 
    fac.reduce = function(coefs) apply(coefs, 2, mean), contr, 
    options = getOption("lsmeans")$lsmeans, weights, trend, ...)
## (used when 'object' is a 'ref.grid' and 'specs' is 'character')
    
\method{lsmeans}{list}(object, specs, ...)
## (used when 'specs' is a 'list')

\method{lsmeans}{formula}(object, specs, contr.list, trend, ...)
## (used when 'specs' is a 'formula')

lstrends(model, specs, var, delta.var = 0.01 * rng, data, 
   transform = c("none", "response"), ...)

lsmobj(bhat, V, levels, linfct, df = NA, post.beta = matrix(NA), ...)

pmmeans(...)
pmtrends(...)
pmmobj(...)
}
\arguments{
\item{object}{
An object of class \code{ref.grid}; or a fitted model object that is supported, such as the result of a call to \code{lm} or \code{lmer}. Many fitted-model objects are supported; see \code{link{models}} for details.}
\item{specs}{
A \code{character} vector specifying the names of the predictors over which LS-means are desired. \code{specs} may also be a \code{formula} or a \code{list} (optionally named) of valid \code{spec}s. Use of formulas is described in the Details section below.}
\item{by}{
A character vector specifying the names of predictors to condition on.}
\item{fac.reduce}{
A function that combines the rows of a matrix into a single vector. This implements the ``marginal averaging'' aspect of least-squares means. The default is the mean of the rows. Typically if it is overridden, it would be some kind of weighted mean of the rows. If \code{fac.reduce} is nonlinear, bizarre results are likely, and LS means will not be interpretable. If the \code{weights} argument is non-missing, \code{fac.reduce} is ignored.}
\item{contr}{
A \code{list} of contrast coefficients to apply to the least-squares means -- or the root name of an \code{.lsmc} function that returns such coefficients. In addition, \code{contr = "cld"} is an alternative way to invoke the \code{\link[lsmeans]{cld}} function. See \code{\link{contrast}} for more details on contrasts. NOTE: \code{contr} is ignored when \code{specs} is a formula.
}
\item{contr.list}{
A named \code{list} of \code{list}s of contrast coefficients, as for \code{contr}. This is used only in the formula method; see Details below. 
}
\item{options}{If non-\code{NULL}, a named \code{list} of arguments to pass to \code{\link{update}}, just after the object is constructed.}
\item{weights}{
Numeric vector, numeric matrix, or character string specifying weights to use in averaging predictions. If a vector, its length must equal the number of predictions to be averaged to obtain each least-squares mean. If a matrix, each row of the matrix is used in turn, wrapping back to the first row as needed.  When in doubt about what is being averaged (or how many), first call with \code{weights = "show.levels"}.)

If a string, it should partially match one of the following:
\describe{
\item{\code{"equal"}}{Use an equally weighted average.}
\item{\code{"proportional"}}{Weight in proportion to the frequencies (in the original data) of the factor combinations that are averaged over.}
\item{\code{"outer"}}{Weight in proportion to each individual factor's marginal frequencies. Thus, the weights for a combination of factors are the outer product of the one-factor margins}
\item{\code{"cells"}}{Weight according to the frequencies of the cells being averaged.}
\item{\code{show.levels}}{This is a convenience feature for understanding what is being averaged over. Instead of a table of LS means, this causes the function to return a table showing the levels that are averaged over, in the order they appear.}
}
Outer weights are like the 'expected' counts in a chi-square test of independence, and will yield the same results as those obtained by proportional averaging with one factor at a time. All except \code{"cells"} uses the same set of weights for each mean. In a model where the predicted values are the cell means, cell weights will yield the raw averages of the data for the factors involved.
Note: If weights were used in fitting the model, then weight totals are used in place of frequencies in these schemes.

If \code{weights} is used, \code{fac.reduce} is ignored. 
}
\item{trend}{
Including this argument is an alternative way of calling \code{lstrends} with \code{trend} as its \code{var} argument and \code{object} as its \code{model}.
}
\item{model}{
A supported model object (\emph{not} a \code{ref.grid}).
}
\item{var}{
Character giving the name of a variable with respect to which a difference quotient of the linear predictors is computed. In order for this to be useful, \code{var} should be a numeric predictor that interacts with at least one factor in \code{specs}. Then instead of computing least-squares means, we compute and compare the slopes of the \code{var} trend over levels of the specified other predictor(s). As in least-squares means, marginal averages are computed when some variables in the reference grid are excluded for the specification.

The user may specify some monotone function of one variable, e.g., \code{var = "log(dose)"}. If so, the chain rule is applied. Note that, in this example, if \code{model} contains \code{log(dose)} as a predictor, we will be comparing the slopes estimated by that model, whereas specifying \code{var = "dose"} would perform a transformation of those slopes.
}
\item{delta.var}{
The value of \emph{h} to use in forming the difference quotient \emph{(f(x+h) - f(x))/h}. Changing it (especially changing its sign) may be necessary to avoid numerical problems such as logs of negative numbers. The default value is 1/100 of the range of \code{var} over the dataset.
}
\item{data}{As in \code{\link{ref.grid}}, you may use this argument to supply the dataset used in fitting the model, for situations where it is not possible to reconstruct the data. Otherwise, leave it missing.}
\item{transform}{In \code{lstrends}, if \code{object} has a response transformation, then specifying \code{transform = "response"} will cause \code{lstrends} to calculate the trends after back-transforming to the response scale. This is done using the chain rule, and standard errors are estimated via the delta method. With \code{transform = "none"} (the default), the trends are calculated on the scale of the linear predictor, without back-transforming it. This argument works similarly to the \code{transform} argument in \code{\link{ref.grid}} (but without a \code{"log"} option), in that the returned object is re-gridded to the new scale (see also \code{\link{regrid}}).}
\item{bhat}{Numeric. Vector of regression coefficients.}
\item{V}{Square matrix. Covariance matrix of \code{bhat}}
\item{levels}{Named list or vector. Levels of factor(s) that define the estimates defined by \code{linfct}. If not a list, we assume one factor named \code{"level"}}
\item{linfct}{Matrix. Linear functions of \code{bhat} for each combination of \code{levels}}
\item{df}{Numeric or function with arguments \code{(x,dfargs)}. If a number, that is used for the degrees of freedom. If a function, it should return the degrees of freedom for \code{sum(x*bhat)}; if additional parameters are needed, include them in \code{\dots} as \code{dfargs} (not abbreviated).}
\item{post.beta}{Matrix whose columns comprise a sample from the posterior distribution of the regression coefficients (so that typically, the column averages will be \code{bhat}). A 1 x 1 matrix of \code{NA} indicates that such a sample is unavailable.}
\item{\dots}{Additional arguments passed to other methods or to \code{\link{ref.grid}}. For example, \code{vcov.} may be used to override the default covariance estimate, and some models allow additional options. Some models require \code{data} to be given explicitly. See the help pages for \code{\link{ref.grid}} and \link{models}. In addition, if the model formula contains references to variables that are not predictors, you must provide a \code{params} argument with a list of their names; see the example below for \code{Oatsq.lm}.}
}


\details{
Least-squares means (also called predicted marginal means) are predictions from a linear model over a \emph{reference grid}, or marginal averages thereof. They have been popularized by \pkg{SAS} (SAS Institute, 2012). The \code{\link{ref.grid}} function identifies/creates the reference grid upon which \code{lsmeans} is based.

For those who prefer the term \dQuote{predicted marginal means}, courtesy wrappers \code{pmmeans}, \code{pmtrends}, and \code{pmmobj} are provided that behave identically to those that start with \code{ls}, except that estimates are relabeled accordingly (e.g., \code{lsmean} becomes \code{pmmean}).

If \code{specs} is a \code{formula}, it should be of the form \code{~ specs}, \code{~ specs | by}, \code{contr ~ specs}, or \code{contr ~ specs | by}. The formula is parsed and the variables therein are used as the arguments \code{specs}, \code{by}, and \code{contr} as indicated. The left-hand side is optional, but if specified it should be the name of a contrast family (e.g., \code{pairwise}) or of a sub-list of \code{contr.list}. Operators like \code{*} or \code{:} are necessary to delineate names in the formulas, but otherwise are ignored.

In the special case where the mean (or weighted mean) of all the predictions is desired, specify \code{specs} as \code{~ 1} or \code{"1"}.

A number of standard contrast families are provided. They can be identified as functions having names ending in \code{.lsmc} -- use 
\preformatted{ls("package:lsmeans", pat=".lsmc")} 
to list them. See the documentation for \code{\link{pairwise.lsmc}} and its siblings for details. You may write your own \code{.lsmc} function for custom contrasts.

The function \code{lsmobj} may be used to construct an object just like one returned by \code{lsmeans} from user-specified coefficients, covariance matrix, levels (or row labels), linear functions for each row, and degrees of freedom. After the object is constructed, it is \code{\link[lsmeans]{update}}ed with any additional arguments in \code{\dots}.
}

\value{
When \code{specs} is a \code{character} vector or one-sided formula, an object of class \code{\link[=lsmobj-class]{lsmobj}}. A number of methods are provided for further analysis, including \code{\link[lsmeans]{summary}}, \code{\link[lsmeans]{confint}}, \code{\link[lsmeans]{test}}, \code{\link[lsmeans]{contrast}}, \code{\link[lsmeans]{pairs}}, and \code{\link[lsmeans]{cld}}.

When \code{specs} is a \code{list} or a \code{formula} having a left-hand side, the return value is an \code{lsm.list} object, which is simply a \code{list} of \code{lsmobj} objects. Methods for \code{lsm.list} objects are the same as those for \code{lsmobj}, but they apply to only one member of the list, determined by its \code{which} argument.

\bold{Side effect:} When \code{object} is a model, a reference grid is constructed and it is saved as \code{.Last.ref.grid} in the user's environment (unless this is disabled via \samp{lsm.option(save.ref.grid = FALSE)}). This makes it possible to check what reference grid was used, or to use it as the \code{object} in future \code{lsmeans} calls (and bypass reconstructing it). Similarly, \code{lstrends} also saves its reference grid (but for predicting difference quotients) as \code{.Last.ref.grid}.
}

\note{If the model formula contains variables that are not predictors (e.g., degree of a polynomial, knots for a spline, etc.), you must add a \code{params} argument to the call}

\note{While using \code{specs} as a two-sided formula or a list is a convenient way to get a lot of results with minimal effort, it can also create confusion when additional arguments are provided, because not all arguments may be applied to all the results produced (see examples). Thus, the safer route is to do things incrementally.}

\note{\code{lsmeans} and its relatives can produce fatal errors or incorrect results with models containing splines (e.g., \code{\link{ns}}) and other smoothers because the required information to reconstruct their basis is not always available. A model with \code{\link{poly}} involving two or more predictors will almost always produce misleading results without any warning; but \code{poly(..., raw = TRUE)} will work correctly.}

\note{For a \code{ref.grid} or \code{lsmobj} object created in \pkg{lsmeans} version 2.10 or earlier, the information needed by the \code{weights} argument is not present; so a message is displayed and averaging is done using \code{fac.reduce}.}


\references{
SAS Institute Inc. (2012)
Online documentation; Shared concepts; LSMEANS statement,
\url{http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_introcom_a0000003362.htm}, accessed August 15, 2012.
}

\author{
Russell V. Lenth
}

\seealso{
\code{\link{ref.grid}}, \code{\link{.Last.ref.grid}}, \code{\link{models}}, \code{\link{pairwise.lsmc}}, \code{\link[multcomp]{glht}}, \code{\link{lsm.options}}
}

\examples{
require(lsmeans)

### Covariance example (from Montgomery Design (8th ed.), p.656)
# Uses supplied dataset 'fiber'
fiber.lm <- lm(strength ~ diameter + machine, data = fiber)

# adjusted means and comparisons, treating machine C as control
( fiber.lsm <- lsmeans (fiber.lm, "machine") )
contrast(fiber.lsm, "trt.vs.ctrlk")
# Or get both at once using
#     lsmeans (fiber.lm, "machine", contr = "trt.vs.ctrlk")


### Factorial experiment
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
( warp.lsm <- lsmeans (warp.lm,  ~ wool | tension,
    options = list(estName = "pred.breaks")) )
pairs(warp.lsm) # remembers 'by' structure
contrast(warp.lsm, method = "poly", by = "wool")


### Unbalanced split-plot example ###
#-- The imbalance is imposed deliberately to illustrate that
#-- the variance estimates become biased
require(nlme)
Oats.lme <- lme(yield ~ factor(nitro) + Variety, 
    random = ~1 | Block/Variety, 
    subset = -c(1,2,3,5,8,13,21,34,55), data = Oats)
(Oats.anal <- lsmeans(Oats.lme, list(poly ~ nitro, pairwise ~ Variety)))


### Issues with lists of specs
test(Oats.anal)                # Uses 1st element by default
confint(Oats.anal, which = 4)  # or confint(Oats.anal[[4]])

# Using 'pmmeans' wrapper ...
pmmeans(warp.lm, ~ wool, 
    options = list(infer = c(TRUE, TRUE), null = 22, side = ">"))

### Weights
# See what's being averaged over in the above
lsmeans(Oats.lme, ~ nitro, cov.reduce = FALSE, weights = "show.levels")

# Give three times the weight to Marvellous
lsmeans(Oats.lme, ~ nitro, cov.reduce = FALSE, weights = c(1,3,1))

# Overall mean
lsmeans(Oats.lme, ~ 1, weights = "equal")
lsmeans(Oats.lme, "1", weights = "cells")


### Model with a quadratic trend for 'nitro'
# Also illustrates use of `params` argument to list non-predictors
deg = 2
Oatsq.lm <- lm(yield ~ Block + poly(nitro, deg) + Variety, data = Oats)
# Predictions at each unique 'nitro' value in the dataset
lsmeans(Oatsq.lm, ~ nitro, cov.reduce = FALSE, params = "deg")


### Trends
fiber.lm <- lm(strength ~ diameter*machine, data=fiber)
# Obtain slopes for each machine ...
( fiber.lst <- lstrends(fiber.lm, "machine", var="diameter") )
# ... and pairwise comparisons thereof
pairs(fiber.lst)

# Suppose we want trends relative to sqrt(diameter)...
lstrends(fiber.lm, ~ machine | diameter, var = "sqrt(diameter)", 
    at = list(diameter = c(20,30)))

# Given summary statistics for 4 cities computed elsewhere, 
# obtain multiple comparisons of their means using the
# Satterthwaite method
ybar <- c(47.6, 53.2, 88.9, 69.8)
s <-    c(12.1, 19.5, 22.8, 13.2)
n <-    c(44,   11,   37,   24)
se2 = s^2 / n
Satt.df <- function(x, dfargs)
    sum(x * dfargs$v)^2 / sum((x * dfargs$v)^2 / (dfargs$n - 1))
city.pmm <- pmmobj(bhat = ybar, V = diag(se2),
    levels = list(city = LETTERS[1:4]), linfct = diag(c(1,1,1,1)),
    df = Satt.df, dfargs = list(v = se2, n = n), estName = "mean")
city.pmm
contrast(city.pmm, "revpairwise")
    
    
# See also many other examples in documentation for 
# 'contrast', 'cld', 'glht', 'lsmip', 'ref.grid', 'MOats',
# 'nutrition', etc., and in the vignettes
}
\keyword{ models }
\keyword{ regression }
\keyword{ htest }
