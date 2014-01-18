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
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object produced by a model-fitting function, e.g. \code{lm}}
  \item{at}{Optional named list of levels for the corresponding variables}
  \item{cov.reduce}{A function, logical value, or a named list of functions, applied to covariates in the underlying data. If a single function, it applies to all covariates. 
  
  If \code{cov.reduce} is a list, then the names of its entries determine which function is used. If no entry is found, the default is \code{mean}. 
  
  If a logical value, \code{TRUE} is equivalent to  \code{cov.reduce=mean}, and \code{FALSE} is equivalent to \code{cov.reduce = function(x) sort(unique(x))}.}
  \item{mult.levs}{A named list of levels for the dimensions of a multivariate response}
}
\details{The reference grid consists of combinations of independent variables over which predictions are made. Least-squares means are defined as these predictions, or marginal averages thereof.
The grid is determined by first reconstructing the data used in fitting the model (see \code{\link{recover.data}}). The default reference grid is determined by the observed levels of any factors, and the results of \code{cov.reduce} for numeric predictors. Any of these may be overridden using \code{at}.

The following model objects are currently supported:
}
\value{An S4 object of class \code{\linkS4Class{ref.grid}}. These objects encapsulate everything needed to do calculations and inferences for least-squares means, and contain nothing that depends on the model-fitting procedure.
}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
c("print.ref.grid", "print.summary.ref.grid", "ref.grid")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
