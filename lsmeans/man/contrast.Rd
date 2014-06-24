\name{contrast}
\alias{contrast}
\alias{contrast.ref.grid}
\alias{contrast.lsm.list}
\alias{test}
\alias{test.ref.grid}
\alias{confint}
\alias{confint.ref.grid}
\alias{pairs}
\alias{pairs.ref.grid}



\title{
Methods for obtaining analyses \code{ref.grid} and \code{lsmobj} objects
}
\description{
These methods provide for analyses of \code{ref.grid} objects, or follow-up analyses of \code{lsmobj} objects: Contrasts, pairwise comparisons, tests, and confidence intervals.
}
\usage{
\method{contrast}{ref.grid}(object, method = "eff", by, adjust, 
    options = getOption("lsmeans")$contrast, ...)
\method{contrast}{lsm.list}(object, ..., which = 1)

\method{test}{ref.grid}(object, parm, ...)

\method{confint}{ref.grid}(object, parm, level = 0.95, ...)

\method{pairs}{ref.grid}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{
An object of class \code{"ref.grid"} or its extension, \code{"lsmobj"}.
}
  \item{method}{
Character value giving the root name of a contast method (e.g. \code{"pairwise"}). Alternatively, a named list of contrast coefficients that must each conform to the number of least-squares means in each \code{by} group. This is just like the \code{contr} argument in \code{\link{lsmeans}}. To identify the available methods, see \preformatted{ls("package:lsmeans", pat=".lsmc")} You may define your own \code{.lsmc} function and use its root name as \code{method}. 
}
  \item{by}{
Character names of variable(s) to be used for ``by'' groups. The contrasts will be evaluated separately for each combination of these variables. If \code{object} was created with by groups, those are used unless overridden. Use \code{by = NULL} to use no by groups at all.
}
  \item{adjust}{
Method to use for adjusting \emph{P} values. This is passed to \code{\link[lsmeans]{summary}}.
}
  \item{options}{If non-\code{NULL}, a named \code{list} of arguments to pass to \code{\link{update}}, just after the object is constructed.}
  \item{parm}{
(This is a required argument for the generic \code{\link[stats]{confint}}) and \code{test} methods, but is ignored by the methods for \code{lsmobj}.)}
  \item{level}{
Numeric value of the desired confidence level.}
  \item{\dots}{
Additional arguments passed to \code{\link[lsmeans]{summary}} or to a contrast function.}
  \item{which}{When \code{object} is a list of \code{lsmobj} objects, this specifies which member of the list is analyzed.}
}
\details{
\code{pairs} is equivalent to \code{contrast} with \code{method = "pairwise"}.

\code{confint} and \code{test} are equivalent to calling \code{\link[lsmeans]{summary}} with \code{infer=c(TRUE,FALSE)} and \code{infer=c(FALSE,TRUE)}, respectively.
}
\value{
All of these methods return an object of class \code{"lsmobj"}, which is an extension of \code{"ref.grid"}. Consequently, they may be used as arguments to other \code{"lsmobj"} or \code{"ref.grid"} methods. The user may, for example, compute contrasts of contrasts, or re-summarize a set of confidence intervals with a different \code{by} grouping or confidence level.

The ``grid'' for the returned value is simply the set of variables that identify the results. For example, \code{contrast}'s return value is a reference grid for one factor named \code{contrast}. 
}
%\references{}
\author{
Russell V. Lenth
}
%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Additional \code{"lsmobj"} methods having their own help pages are \code{\link[lsmeans]{cld}} and \code{\link[lsmeans]{glht}}. Also, the \code{\link[lsmeans]{summary}} and other methods for \code{"ref.grid"} objects also work for \code{"lsmobj"} objects.
}
\examples{
require(lsmeans)
warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
warp.lsm <- lsmeans(warp.lm, ~ tension | wool)

# Polynomial contrasts of tension, by wool
(warp.pl <- contrast(warp.lsm, "poly"))
# Same results with a different adjustment
summary(warp.pl, adjust = "fdr")

# Compare the two contrasts for each degree
contrast(warp.pl, "revpairwise", by = "contrast")

# User-provided contrasts, ignoring the previous by grouping
contrast(warp.lsm, 
    list(c1=c(1,0,0,-1,0,0), c2=c(1,1,1,-1,-1,-1)/3), 
    by = NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
\keyword{ htest }

