\name{glht}
\alias{glht}
\alias{lsm}
\alias{glht.lsmobj}
\alias{as.glht}
\alias{as.glht.lsmobj}
\alias{summary.glht.list}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\pkg{lsmeans} support for \code{glht}
}
\description{
These functions and methods provide an interface between \pkg{lsmeans} and the \code{\link[multcomp]{glht}} function for simultaneous inference in the \pkg{multcomp} package.
}
\usage{
\method{as.glht}{lsmobj}(object, ...)

\method{summary}{glht.list}(object, ...)

lsm(...)
}
\arguments{
  \item{object}{
An object of the required class.
}
  \item{\dots}{
Additional arguuments to other methods.
}
}
\details{
\code{lsm} is meant to be called only \emph{from} \code{"glht"} as its second (\code{linfct}) argument. It works similarly to \code{\link[multcomp]{mcp}} except with \code{specs} (and optionally \code{by} and \code{contr} arguments) provided as in a call to \code{\link{lsmeans}}. 

When there is a non-\code{NULL} \code{by} variable (either explicitly or implicitly), each ``by'' group is passed separately to \code{glht} and returned as a \code{list} of \code{"glht"} objects. For convenience, this is classed as \code{"glht.list"} and a \code{summary} method is provided.
}
\value{
\code{as.glht} returns an object of class \code{\link[multcomp]{glht}}, or of class \code{glht.list} if \code{by} is non-\code{NULL}. The latter is simply a list of \code{glht} objects but it has its own \code{summary} method which returns a list of \code{summary.glht} objects. 
}
%\references{}
\author{
Russell V. Lenth
}
\note{
There is also a \code{glht} method for class \code{lsmobj}, but it is far preferable to use \code{as.glht} instead, as \code{model} is redundant. 
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
