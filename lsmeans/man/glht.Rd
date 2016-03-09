\name{glht}
\alias{glht}
\alias{lsm}
\alias{glht.ref.grid}
\alias{as.glht}
\alias{as.glht.ref.grid}
\alias{as.glht.lsm.list}
\alias{coef.glht.list}
\alias{confint.glht.list}
\alias{plot.glht.list}
\alias{summary.glht.list}
\alias{vcov.glht.list}
\alias{pmm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\pkg{lsmeans} support for \code{glht}
}
\description{
These functions and methods provide an interface between \pkg{lsmeans} and the \code{\link[multcomp]{glht}} function for simultaneous inference in the \pkg{multcomp} package.
}
\usage{
\method{as.glht}{ref.grid}(object, ...)
\method{as.glht}{lsm.list}(object, ..., which = 1)

\method{coef}{glht.list}(object, ...)
\method{confint}{glht.list}(object, ...)
\method{plot}{glht.list}(x, ...)
\method{summary}{glht.list}(object, ...)
\method{vcov}{glht.list}(object, ...)

lsm(...)
pmm(...)
}
\arguments{
  \item{object, x}{
An object of the required class.
}
  \item{which}{Numeric index of which element of the \code{lsm.list} to use.}
  \item{\dots}{
Additional arguuments to other methods.
}
}
\details{
\code{lsm} (and \code{pmm}, which is identical) are meant to be called only \emph{from} \code{"glht"} as its second (\code{linfct}) argument. It works similarly to \code{\link[multcomp]{mcp}} except with \code{specs} (and optionally \code{by} and \code{contr} arguments) provided as in a call to \code{\link{lsmeans}} or \code{pmmeans}. 

When there is a non-\code{NULL} \code{by} variable (either explicitly or implicitly), each ``by'' group is passed separately to \code{glht} and returned as a \code{list} of \code{"glht"} objects. For convenience, this is classed as \code{"glht.list"}, and appropriate methods \code{coef}, \code{confint}, \code{plot}, \code{summary}, and \code{vcov} are provided.
}
\value{
\code{as.glht} returns an object of class \code{\link[multcomp]{glht}}, or of class \code{glht.list} if \code{by} is non-\code{NULL}. The latter is simply a list of \code{glht} objects, and the provided methods \code{coef}, \code{confint}, \code{plot}, \code{summary}, and \code{vcov} simply \code{\link{lapply}} the corresponding methods for class \code{"glht"}. 
}
%\references{}
\author{
Russell V. Lenth
}
\note{
There is also a \code{glht} method for class \code{ref.grid}, but it is far preferable to use \code{as.glht} instead, as its \code{model} argument is redundant. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lsmeans}}, \code{\link[multcomp]{glht}}
}
\examples{
require(lsmeans)
require(multcomp)

warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)

# Using 'lsm'
summary(glht(warp.lm, lsm(pairwise ~ tension | wool)))

# Same, but using an existing 'lsmeans' result
warp.lsmobj <- lsmeans(warp.lm, ~ tension | wool)
summary(as.glht(pairs(warp.lsmobj)))

# Same contrasts, but treat as one family
summary(as.glht(pairs(warp.lsmobj), by = NULL))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ models }
\keyword{ htest }
