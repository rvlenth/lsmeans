\name{lsmeans}
\alias{lsmeans}
\alias{lsmeans.formula}
\alias{lsmeans.list}
\alias{lsmeans.character}
\alias{lsmeans.character.ref.grid}
\alias{lstrends}

\title{Least-squares means}
\description{
Compute least-squares means for specified factors or factor combinations in a linear model,
and optionally comparisons or contrasts among them.
}
\usage{
\method{lsmeans}{character}(object, specs, ...)
## (used when 'specs' is 'character')

\method{lsmeans}{character.ref.grid}(object, specs, by = NULL, 
    fac.reduce = function(coefs) apply(coefs, 2, mean), contr, 
    options = getOption("lsmeans")$lsmeans, ...)
## (used when 'object' is a 'ref.grid' and 'specs' is 'character')
    
\method{lsmeans}{list}(object, specs, ...)
## (used when 'specs' is a 'list')

\method{lsmeans}{formula}(object, specs, contr.list, trend, ...)
## (used when 'specs' is a 'formula')

lstrends(model, specs, var, delta.var = 0.01 * rng, data, ...)
}
\arguments{
\item{object}{
An object of class \code{ref.grid}; or a fitted model object that is supported, such as the result of a call to \code{lm} or \code{lmer}.}
\item{specs}{
A \code{character} vector specifying the names of the predictors over which LS-means are desired. \code{specs} may also be a \code{formula} or a \code{list} (optionally named) of valid \code{spec}s. Use of formulas is described in the Details section below.}
\item{by}{
A character vector specifying the names of predictors to condition on.}
\item{fac.reduce}{
A function that combines the rows of a matrix into a single vector. This implements the ``marginal averaging'' aspect of least-squares means. The default is the mean of the rows. Typically if it is overridden, it would be some kind of weighted mean of the rows. If \code{fac.reduce} is nonlinear, bizarre results are likely, and LS means will not be interpretable.}
\item{contr}{
A \code{list} of contrast coefficients to apply to the least-squares means -- or the root name of an \code{.lsmc} function that returns such coefficients. In addition, \code{contr = "cld"} is an alternative way to invoke the \code{\link[lsmeans]{cld}} function. See \code{\link{contrast}} for more details on contrasts.
}
\item{contr.list}{
A named \code{list} of \code{list}s of contrast coefficients, as for \code{contr}. This is used only in the formula method; see Details below. 
}
\item{options}{If non-\code{NULL}, a named \code{list} of arguments to pass to \code{\link{update}}, just after the object is constructed.}
\item{trend}{
Including this argument is an alternative way of calling \code{lstrends} with it as its \code{var} argument.
}
\item{model}{
A supported model object.
}
\item{var}{
Character giving the name of a variable with respect to which a difference quotient of the linear predictors is computed. In order for this to be useful, \code{var} should be a numeric predictor that interacts with at least one factor in \code{specs}. Then instead of computing least-squares means, we compute and compare the slopes of the \code{var} trend over levels of the specified other predictor(s). As in least-squares means, marginal averages are computed when some variables in the reference grid are excluded for the specification.

The user may specify some monotone function of one variable, e.g., \code{var = "log(dose)"}. If so, the chain rule is applied. Note that, in this example, if \code{model} contains \code{log(dose)} as a predictor, we will be comparing the slopes estimated by that model, whereas specifying \code{var = "dose"} would perform a transformation of those slopes.
}
\item{delta.var}{
The value of \emph{h} to use in forming the difference quotient \emph{(f(x+h) - f(x))/h}. Changing it (especially changing its sign) may be necessary to avoid numerical problems such as logs of negative numbers. The default value is 1/100 of the range of \code{var} over the dataset.
}
\item{data}{As in \code{\link{ref.grid}}, you may use this argument to supply the dataset used in fitting the model, for situations where it is not possible to reconstruct the data. Otherwise, leave it missing.}
\item{\dots}{Additional arguments passed to other methods or to \code{\link{ref.grid}}.}
}


\details{
Least-squares means are predictions from a linear model over a \emph{reference grid}, or marginal averages thereof. They have been popularized by \pkg{SAS} (SAS Institute, 2012). The \code{\link{ref.grid}} function identifies/creates the reference grid upon which \code{lsmeans} is based.

If \code{specs} is a \code{formula}, it should be of the form \code{contr ~ specs | by}. The formula is parsed and then used as the arguments \code{contr}, \code{specs}, and \code{by} as indicated. The left-hand side is optional, but if specified it should be the name of a contrast family (e.g., \code{pairwise}) or of a sub-list of \code{contr.list}. Operators like \code{*} or \code{:} are necessary to delineate names in the formulas, but otherwise are ignored.

A number of standard contrast families are provided. They can be identified as functions having names ending in \code{.lsmc} -- use \preformatted{ls("package:lsmeans", pat=".lsmc")} to list them. See the documentation for \code{\link{pairwise.lsmc}} and its siblings for details. You may write your own \code{.lsmc} function for custom contrasts.

For models fitted using the \pkg{lme4} package, degrees of freedom are obtained using the Kenward-Roger (1997) method as implemented in the package \pkg{pbkrtest}, if it is installed. If \pkg{pbkrtest} is not installed, the degrees of freedom are set to \code{NA} and asymptotic results are displayed.
}

\value{
An object of class \code{\link[=lsmobj-class]{lsmobj}} -- except when \code{specs} is a \code{list} or a \code{formula} having a left-hand side, a \code{list} of \code{slmobj} objects. A number of methods are provided for further analysis, including \code{\link[lsmeans]{summary}}, \code{\link[lsmeans]{confint}}, \code{\link[lsmeans]{test}}, \code{\link[lsmeans]{contrast}}, \code{\link[lsmeans]{pairs}}, and \code{\link[lsmeans]{cld}}.
}

%\note{}

\references{
Kenward, M.G. and Roger, J.H. (1997)
Small sample inference for fixed effects from restricted maximum likelihood,
\emph{Biometrics}, 53, 983--997.

SAS Institute Inc. (2012)
Online documentation; Shared concepts; LSMEANS statement,
\url{http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_introcom_a0000003362.htm}, accessed August 15, 2012.
}

\author{
Russell V. Lenth
}

\seealso{
\code{\link{pairwise.lsmc}}, \code{\link[multcomp]{glht}}
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
lsmeans(Oats.lme, list(poly ~ nitro, pairwise ~ Variety))

# Model with a quadratic trend for 'nitro'
Oatsq.lme <- update(Oats.lme, . ~ nitro + I(nitro^2) + Variety)
# Predictions at each unique 'nitro' value in the dataset
lsmeans(Oatsq.lme, ~ nitro, cov.reduce = FALSE)

# Trends
fiber.lm <- lm(strength ~ diameter*machine, data=fiber)
# Obtain slopes for each machine ...
( fiber.lst <- lstrends(fiber.lm, "machine", var="diameter") )
# ... and pairwise comparisons thereof
pairs(fiber.lst)

# Suppose we want trends relative to sqrt(diameter)...
lstrends(fiber.lm, ~ machine | diameter, var = "sqrt(diameter)", 
    at = list(diameter = c(20,30)))
    
# See also many other examples in documentation for 
# 'contrast', 'cld', 'glht', 'lsmip', 'ref.grid', 'MOats',
# 'nutrition', etc., and in the vignettes
}
\keyword{ models }
\keyword{ regression }
\keyword{ htest }
