\name{summary}
\alias{summary}
\alias{summary.ref.grid}
\alias{str.ref.grid}
\alias{print.ref.grid}
\alias{print.summary.ref.grid}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for \code{ref.grid} objects
}
\description{
Use these methods to summarize, print, or examine objects of class \code{"ref.grid"}. They also apply to the class \code{"lsmobj"}, which is an extension of \code{"ref.grid"}.
}
\usage{
\method{summary}{ref.grid}(object, infer, level, adjust, by, ...)

\method{str}{ref.grid}(object, ...)

\method{print}{ref.grid}(x, ...)
\method{print}{summary.ref.grid}(x, ..., digits = NULL, quote = FALSE, right = TRUE)
}

\arguments{
  \item{object}{
An object of class \code{"ref.grid"}.
}
  \item{infer}{
A vector of two logical values. The first determines whether confidence intervals are displayed, and the second determines whether \emph{t} tests and \emph{P} values are displayed. If only one value is provided, it is used for both. 
}
  \item{level}{
Confidence level for confidence intervals, if \code{infer[1]} is \code{TRUE}.
}
  \item{adjust}{
Character value naming the method used to adjust \emph{P} values, if \code{infer[2]} is \code{TRUE}. 
\code{"tukey"} computes p values using the Studentized range distribution with the number of means in the family;
\code{"sidak"} replaces each p value by \code{1 - (1 - p)^c}, where c is the number of contrasts;
\code{"scheffe"} computes p values from the F distribution, according to the Scheffe criterion;
\code{"none"} makes no adjustments to the p values.
In some cases, these adjustments are only approximate, especially when the degrees of freedom or standard errors vary
greatly within the family of tests. For more accurate adjustments, use the \code{\link[lsmeans]{glht}} method for \pkg{multcomp}.

In addition, the Bonferroni-inequality-based adjustment methods in \code{\link{p.adjust}} are available.
}
  \item{by}{
Character name(s) of variables to use for grouping. This affects the family of tests considered in adjusted \emph{P} values. The printed display of the summary is grouped by the \code{by} variables.
}
  \item{x}{
The object to be printed.
}
  \item{\dots, digits, quote, right}{Additional arguments passed to \code{\link{print.data.frame}}.}
}
\details{
The \code{misc} slot in \code{object} contains default values for these parameters. The defaults vary depending on the code that created the object.
}
\value{
The \code{summary} method for \code{"ref.grid"} objects returns An object of class \code{"summary.ref.grid"}, which extends \code{"data.frame"}.
}
%%\references{}
\author{
Russell V. Lenth
}
%%%\note{}

\seealso{
Methods for the closely related \code{"lsmobj"} class can be found in \code{\link[lsmeans]{contrast}}
}
\examples{
require(lsmeans)
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.rg <- ref.grid(warp.lm)
str(warp.rg)

summary(warp.rg)

summary(warp.rg, by = "wool", 
        infer = c(TRUE, FALSE), level = .90)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
