\name{make.tran}
\alias{make.tran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Response transformations
}
\describe{
Create the needed information to perform transformations of the response variable, including inverting the transformation and estimating variances of back-transformed predictions via the delta method. \code{make.tran} is similar to \code{\link{make.link}}, but it covers additional transformations. The result can be used as an environment in which the model is fitted, or as the \code{tran} argument in \code{\link{update.ref.grid}} (when the given transformation was already applied in an existing model).
}
\usage{
make.tran(type = c("genlog", "power", "boxcox", "arcsine"), param = 1)

# Additional named transformations not requiring make.tran:
#   "reciprocal"  "log2"   "log10"  "arcsine"
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
The name of the transformation. See Details.
}
  \item{param}{
A parameter for the transformation. See Details.
}
}
\details{
The functions \code{\link{lsmeans}}, \code{\link{ref.grid}}, and related ones automatically detect response transformations that are recognized by examining the model formula. Certain standard ones such as \code{log} and \code{sqrt} -- especially ones commonly also used as link functions in generalized linear models --  are recognized by name and supported via \code{\link{make.link}}. The \pkg{lsmeans} package recognizes a few others (such as \code{log10}) by name, as listed in Usage. 

Support of named transformations is not enough for more general transformations; so a \code{ref.grid} or \code{lsmobj} object may also be modified via \code{\link{update}} so that its \code{tran} component is a list of functions like those returned by \code{make.link}. The \code{make.tran} function extends the capabilities of \code{make.link}, and allows for a wider selection of response transformations. If the returned list, say \cdoe{tran.result}, is used as an environment within which the model is fitted, and the response is transformed using the function \code{linkfun}, \code{ref.grid} recognizes this and automatically substitutes \code{tran.result} as its \code{tran} component.

It is also possible that an existing fitted model has a response transformation available as a supported named transformation or supported by \code{make.tran}. In that case, it is not necessary to re-fit the model; instead, update the result of \code{ref.grid} or \code{lsmeans} with the appropriate transformation name or list.

Most of the transformations available in "make.tran" require a parameter, specified in \code{param}; we use \eqn{p} to denote this parameter, and \eqn{y} to denote the response variable, in subsequent expressions. The \code{type} argument specifies the following transformations:
\describe{
\item{\code{genlog}}{Generalized logarithmic transformation: \eqn{log(y + p)} }
\item{\code{power}}{Power transformation: \eqn{y^p}}
\item{\code{boxcox}}{The Box-Cox transformation (unscaled by the geometric mean): \eqn{(y^p - 1) / p} when \eqn{p} is nonzero, \eqn{log(y)} when \eqn{p=0}}
\item{\code{arcsine}}{Arcsin-square-root transformation: \eqn{sin^(-1)y^{1/2)}. This requires no parameters, and is also available as a named transformation.}
}
}
\value{
A \code{list} having at least the same elements as that returned by \code{\link{make.link}}. The \code{linkfun} component is the transformation itself. 
}
%%%\references{}
\author{
Russell V. Lenth
}
\note{
We modify certain \code{make.link} results in transformations where there is a restriction on valid prediction values, so that reasonable inverse predictions are obtained for no matter what. For example, if a \code{sqrt} transformation was used but a predicted value is negative, the inverse transformation is zero rather than the square of the prediction. A side effect of this is that it is possible for one or both confidence limits, or even a standard error, to be zero.
}

\seealso{
  \code{\link{make.link}}, \code{\link[lsmeans]{update}}
}
\examples{
require("lsmeans")

# Fit a model using an oddball transformation:
bctran <- make.tran("boxcox", 0.368)
warp.bc <- with(bctran, 
    lm(linkfun(breaks) ~ wool * tension, data = warpbreaks))
# Obtain back-transformed LS means:    
lsmeans(warp.bc, type = "response")

\dontrun{
# An existing model 'mod' was fitted with a log(y + 1) transformation...
mod.rg <- update(ref.grid(mod), tran = make.tran("genlog", 1))
lsmeans(mod.rg, "treatment")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
