\name{glht}
\alias{glht}
\alias{lsm}
\alias{glht.lsmobj}
\alias{summary.glht.list}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\pkg{lsmeans} support for \code{glht}
}
\description{
These functions and methods provide an interface between \pkg{lsmeans} and the \code{\link[multcomp]{glht}} function for simultaneous inference in the \pkg{multcomp} package.
}
\usage{
\method{glht}{lsmobj}(model, linfct, by, ...)

\method{summary}{glht.list}(object, ...)

lsm(...)
}
\arguments{
  \item{model}{
A model object supported by \code{glht}
}
  \item{linfct}{
An object of class \code{\link[=lsmobj-class]{lsmobj}} or \code{lsmlf}
}
  \item{by}{
Character name(s) of grouping variables; or \code{NULL} to disable ``by'' variables already there.
}
  \item{object}{
An object of class \code{"glht.list"} (which is simply a list of \code{"glht"} objects).
}
  \item{\dots}{
Additional arguuments to other methods.
}
}
\details{
\code{lsm} is meant to be called only \emph{from} \code{"glht"} as its second (\code{linfct}) argument. It works similarly to \code{\link[multcomp]{mcp}} except with \code{specs}, (and optionally \code{by} and \code{contr} arguments) as in a call to \code{\link{lsmeans}}. 

When there is a non-\code{NULL} \code{by} variable (either explicitly or implicitly), each ``by'' group is passed separately to \code{glht} and returned as a \code{list} of \code{"glht"} objects. For convenience, this is classed as \code{"glht.list"} and a \code{summary} method is provided.
}
\value{
\code{glht} returns an object of class \code{\link[multcomp]{glht}}, or if \code{by} is non-\code{NULL}, a list of \code{glht} objects. 
}
%\references{}
\author{
Russell V. Lenth
}
%\note{}

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

# Same, passing an 'lsmobj' object
warp.lsmobj <- lsmeans(warp.lm, ~ tension | wool)
summary(glht(warp.lm, pairs(warp.lsmobj)))

# Same contrasts, but treat as one family
summary(glht(warp.lm, pairs(warp.lsmobj), by = NULL))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
\keyword{ models }
\keyword{ htest }
