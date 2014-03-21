\name{summary}
\alias{summary}
\alias{summary.ref.grid}
\alias{str.ref.grid}
\alias{predict.ref.grid}
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
\method{summary}{ref.grid}(object, infer, level, adjust, by, 
    type = c("link","response","lp","linear"), ...)

\method{predict}{ref.grid}(object, type = c("link","response","lp","linear"), ...)

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
Character value naming the method used to adjust \eqn{p} values or confidence limits. See Details.
}
  \item{by}{
Character name(s) of variables to use for grouping. This affects the family of tests considered in adjusted \emph{P} values. The printed display of the summary is grouped by the \code{by} variables.
}
  \item{type}{
Type of prediction desired. This only has an effect if there is a known transformation or link function. \code{"response"} specifies that the inverse transformation be applied. All others are equivalent, and request that results be shown for the linear predictor.
}
  \item{x}{
The object to be printed.
}
  \item{\dots, digits, quote, right}{Additional arguments passed to \code{\link{print.data.frame}}.}
}
\details{
The \code{misc} slot in \code{object} contains default values for \code{by}, \code{infer}, \code{level}, and \code{adjust}. These defaults vary depending on the code that created the object.

With \code{type="response"}, the transformation assumed can be found in \samp{object@misc$tran}, and its label, for the summary is in \samp{object@misc$inv.lbl}. At this time, \code{tran} must be one of the named transformations valid for \code{\link{make.link}}. Any \eqn{t} or \eqn{z} tests are still performed on the scale of the linear predictor, not the inverse-transformed one. Similarly, confidence intervals are computed on the linear-predictor scale, then inverse-transformed.

As applied to \eqn{p} values, \code{adjust} has the following effects:
\code{"tukey"} computes \eqn{p} values using the Studentized range distribution with the number of means in the family;
\code{"sidak"} replaces each \eqn{p} value by \eqn{1 - (1 - p)^c}, where \eqn{c} is the number of contrasts;
\code{"scheffe"} computes \eqn{p} values from the \eqn{F} distribution, according to the Scheffe criterion;
\code{"none"} makes no adjustments to the \eqn{p} values.
In addition, the Bonferroni-inequality-based adjustment methods in \code{\link{p.adjust}} are available.

Adjustments for confidence limits are available only for \code{adjust} equal to \code{"tukey"}, \code{"sidak"}, \code{"bonferroni"}, and \code{"none"}. For any other \code{adjust}, confidence limits are left unadjusted and a message is displayed to that effect.

In some cases, confidence and \eqn{p}-value adjustments are only approximate -- especially when the degrees of freedom or standard errors vary
greatly within the family of tests. For more accurate adjustments, use the \code{\link[lsmeans]{glht}} method for \pkg{multcomp}.
}
\value{
The \code{summary} method for \code{"ref.grid"} objects returns an object of class \code{"summary.ref.grid"}, which extends \code{"data.frame"}.
}
%%\references{}
\author{
Russell V. Lenth
}
%%%\note{}

\seealso{
Methods for the closely related \code{"lsmobj"} class can be found in \code{\link[lsmeans]{contrast}}, \code{\link[lsmeans]{cld}}, and \code{\link[lsmeans]{glht}}.
}
\examples{
require(lsmeans)
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.rg <- ref.grid(warp.lm)
str(warp.rg)

summary(warp.rg)

summary(warp.rg, by = "wool", 
        infer = c(TRUE, FALSE), level = .90, adjust = "sidak")

# Transformed response
sqwarp.rg <- ref.grid(update(warp.lm, sqrt(breaks) ~ .))
summary(sqwarp.rg)

# Back-transformed results - compare with summary of 'warp.rg'
summary(sqwarp.rg, type = "response")

# Logistic regression
# Reshape the Titanic data
Titan <- do.call("expand.grid", dimnames(Titanic)[-4])
Titan$Died <- matrix(Titanic, ncol=2)
Titan.glm <- glm(Died ~ (Class + Sex + Age)^2, 
    family = binomial, data = Titan)
Titan.lsm <- lsmeans(Titan.glm, ~ Class|Sex, at = list(Age="Adult"))
summary(Titan.lsm, type="response")
summary(pairs(Titan.lsm), type="response")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
