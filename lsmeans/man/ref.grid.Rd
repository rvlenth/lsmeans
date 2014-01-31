\name{ref.grid}
\alias{ref.grid}


% handy for extra documentation...
\alias{ref.grid.object} 
\alias{summary.ref.grid.object}

\title{
Create a reference grid from a fitted model
}
\description{
Using a fitted model object, determine a reference grid for which least-squares means are defined. The resulting \code{ref.grid} object encapsulates all the information needed to calculate LS means and make inferences on them.
}
\usage{
ref.grid(object, at, cov.reduce = mean, mult.levs)
}

\arguments{
  \item{object}{An object produced by a supported model-fitting function, such as \code{lm}. See \code{methods("lsm.basis")} to get an idea of what model objects are supported.}
  \item{at}{Optional named list of levels for the corresponding variables}
  \item{cov.reduce}{A function, logical value, or a named list of functions, applied to covariates in the underlying data. If a single function, it applies to all covariates. 
  
  If \code{cov.reduce} is a list, then the names of its entries determine which function is used. If no entry is found, the default is \code{mean}. 
  
  If a logical value, \code{TRUE} is equivalent to  \code{cov.reduce=mean}, and \code{FALSE} is equivalent to \code{cov.reduce = function(x) sort(unique(x))}.
} % end of \item{cov.reduce}
  \item{mult.levs}{A named list of levels for the dimensions of a multivariate response. If there is more than one element, the combinations of levels are used, in \code{\link{expand.grid}} order. The (total) number of levels must match the number of dimensions. If there is a multivariate response and \code{mult.levs} is missing, the default name is \code{rep.meas} with levels \code{1,2,...}}
} %--- end of \arguments

\details{The reference grid consists of combinations of independent variables over which predictions are made. Least-squares means are defined as these predictions, or marginal averages thereof.
The grid is determined by first reconstructing the data used in fitting the model (see \code{\link{recover.data}}). The default reference grid is determined by the observed levels of any factors, and the results of \code{cov.reduce} for numeric predictors. Any of these may be overridden using \code{at}.

Ability to support a particular class of \code{object} depends on the existence of \code{recover.data} and \code{lsm.basis} methods -- see \link{extending-lsmeans} for details. The call \code{methods("recover.data")} will help identify these.
} % --- end of details

\value{An S4 object of class \code{"ref.grid"} (see \code{\link{ref.grid.class}}). These objects encapsulate everything needed to do calculations and inferences for least-squares means, and contain nothing that depends on the model-fitting procedure.
}


%%%\references{}

\author{
Russell V. Lenth
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link[lsmeans]{summary}} and other methods for the returned objcts. Reference grids are fundamental to \code{\link{lsmeans}}. Click here for more on the \code{\link[=ref.grid-class]{ref.grid}} class. 
}

\examples{
require(lsmeans)

fiber.lm <- lm(strength ~ machine*diameter, data = fiber)
ref.grid(fiber.lm)
summary(ref.grid(fiber.lm, at = list(diameter = c(15,25))))

# Multivariate example
MOats.lm = lm(yield ~ Block + Variety, data = MOats)
ref.grid(MOats.lm, mult.levs = list(nitro = c(0,.2,.4,.6)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ models }% __ONLY ONE__ keyword per line
