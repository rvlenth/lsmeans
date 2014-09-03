\name{summary}
\alias{summary}
\alias{summary.ref.grid}
\alias{str.ref.grid}
\alias{predict.ref.grid}
\alias{print.ref.grid}
\alias{print.summary.ref.grid}
\alias{plot.lsmobj}
\alias{plot.summary.ref.grid}
\alias{vcov.ref.grid}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for \code{ref.grid} objects
}
\description{
Use these methods to summarize, print, plot, or examine objects of class \code{"ref.grid"}. They also apply to the class \code{"lsmobj"}, which is an extension of \code{"ref.grid"}.
}
\usage{
\method{summary}{ref.grid}(object, infer, level, adjust, by, type, df, ...)

\method{predict}{ref.grid}(object, type, ...)

\method{str}{ref.grid}(object, ...)

\method{print}{ref.grid}(x, ...)
\method{print}{summary.ref.grid}(x, ..., digits = NULL, quote = FALSE, right = TRUE)

\method{plot}{lsmobj}(x, y, type, intervals = TRUE, comparisons = FALSE,
    alpha = 0.05, adjust = "tukey", int.adjust = "none", ...)
\method{plot}{summary.ref.grid}(x, y, horizontal = TRUE, 
    xlab, ylab, layout, ...)

\method{vcov}{ref.grid}(object, ...)
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
Character value naming the method used to adjust \eqn{p} values or confidence limits; or to adjust comparison arrows in \code{plot}. See Details.
}
  \item{by}{
Character name(s) of variables to use for grouping. This affects the family of tests considered in adjusted \emph{P} values. The printed display of the summary is grouped by the \code{by} variables.
}
  \item{type}{
Type of prediction desired. This only has an effect if there is a known transformation or link function. \code{"response"} specifies that the inverse transformation be applied. Other valid values are \code{"link"}, \code{"lp"}, and \code{"linear"}; these are equivalent, and request that results be shown for the linear predictor. The default is \code{"link"}, unless the \code{"predict.type"} option is in force; see \code{\link{lsm.options}}.
}
  \item{df}{
If non-missing a constant number of degrees of freedom to use in constructing confidence intervals and \emph{P} values (\code{NA} specifies asymptotic results).
}
  \item{x}{
The object to be printed or plotted.
}
  \item{y}{This argument is ignored.}
  \item{horizontal}{Determines orientatiion of plotted confidence intervals.}
  \item{intervals}{If \code{TRUE}, confidence intervals are plotted for each estimate}
  \item{comparisons}{If \code{TRUE}, \dQuote{comparison arrows} are added to the plot, in such a way that the degree to which arrows overlap reflects as much as possible the significance of the comparison of the two estimates.}
  \item{alpha, int.adjust}{The \code{alpha} argument to use in constructing comparison arrows. \code{int.adjust} may be used to set the \code{adjust} argument for the confidence intervals (use \code{adjust} to set the adjust method for the comparison arrows).}
  \item{\dots, digits, quote, right, xlab, ylab, layout}{For summaries, these are additional arguments passed to other methods including \code{\link{print.data.frame}}, \code{\link{update}}, or \code{\link{dotplot}} as appropriate. If not specified, appropriate defaults are used. For example, the default \code{layout} is one column of horizontal panels or one row of vertical panels.}
}
\details{
The \code{misc} slot in \code{object} contains default values for \code{by}, \code{infer}, \code{level}, \code{adjust}, and \code{type}. These defaults vary depending on the code that created the object. The \code{\link{update}} method may be used to change these defaults. In addition, any options set using \samp{lsm.options(summary=...)} will trump those stored in the object's \code{misc} slot.

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

When the model is rank-deficient, each row \code{x} of \code{object}'s \code{linfct} slot is each checked for estimability. If \code{sum(x*bhat)} is found to be non-estimable, then an \code{NA} is displayed for the estimate (as well as any associated statistics). This check is performed using the orthonormal basis \code{N} in the \code{nbasis} slot for the null space of the rows of the model matrix. Estimability fails when \eqn{||Nx||^2 / ||x||^2} exceeds \code{tol}, which by default is \code{1e-8}. You may change it via \code{\link{lsm.options}} by setting \code{estble.tol} to the desired value.

The \code{plot} method for \code{"lsmobj"} or \code{"summary.ref.grid"} objects (but not \code{"ref.grid"} objects themselves) produces a plot displaying confidence intervals for the estimates. If any \code{by} variables are in force, the plot is divided into separate panels. These functions use the \code{\link[lattice]{dotplot}} function, and thus require that the \pkg{lattice} package be installed. For \code{"summary.ref.grid"} objects, the \code{\dots} arguments in \code{plot} are passed \emph{only} to \code{dotplot}, whereas for \code{"lsmobj"} objects, the object is updated using \code{\dots} before summarizing and plotting.

In plots with \code{comparisons = TRUE}, the resulting arrows are only approximate, and in some cases may fail to accurately reflect the pairwise comparisons of the estimates -- especially when estimates having large and small standard errors are intermingled in just the wrong way.
}

\value{
The \code{summary} method for \code{"ref.grid"} objects returns an object of class \code{"summary.ref.grid"}, which extends \code{"data.frame"}. \code{plot} returns an object of class \code{"trellis"}. \code{vcov} returns the covariance matrix of the object's \code{linfct} slot.
}
%%\references{}
\author{
Russell V. Lenth
}
%\note{}

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

# Plot 90% CIs on the response scale
plot(Titan.lsm, type = "response", level = .90, 
     xlab = "Predicted probability of drowning")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
