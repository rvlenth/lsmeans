%    Copyright (c) 2012-2016 Russell V. Lenth                                %

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
\alias{regrid}
\alias{rbind}
\alias{rbind.ref.grid}
\alias{as.mcmc}
\alias{as.mcmc.ref.grid}
\alias{[.ref.grid}
\alias{xtable.ref.grid}
\alias{xtable.summary.ref.grid}
\alias{print.xtable.lsm}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for \code{ref.grid} objects
}
\description{
Use these methods to summarize, print, plot, or examine objects of class \code{"ref.grid"}. They also apply to the class \code{"lsmobj"}, which is an extension of \code{"ref.grid"}.
}
\usage{
\method{summary}{ref.grid}(object, infer, level, adjust, by, type, df, 
    null, delta, side, ...)

\method{predict}{ref.grid}(object, type, ...)

\method{str}{ref.grid}(object, ...)

\method{rbind}{ref.grid}(..., deparse.level = 1, adjust = "mvt")

\method{[}{ref.grid}(x, i, adjust = "mvt", ...)

\method{print}{ref.grid}(x, ...)
\method{print}{summary.ref.grid}(x, ..., digits = NULL, quote = FALSE, right = TRUE)

\method{xtable}{ref.grid}(x, caption = NULL, label = NULL, align = NULL, 
    digits = 4, display = NULL, auto = FALSE, ...)
\method{xtable}{summary.ref.grid}(x, caption = NULL, label = NULL, align = NULL, 
    digits = 4, display = NULL, auto = FALSE, ...)
\method{print}{xtable.lsm}(x, type = getOption("xtable.type", "latex"),
    include.rownames = FALSE, sanitize.message.function = footnotesize,
    ...)

\method{plot}{lsmobj}(x, y, type, intervals = TRUE, comparisons = FALSE,
    alpha = 0.05, adjust = "tukey", int.adjust, ...)
\method{plot}{summary.ref.grid}(x, y, horizontal = TRUE, 
    xlab, ylab, layout, ...)

\method{vcov}{ref.grid}(object, ...)

regrid (object, transform = c("response", "log", "none"), 
    inv.log.lbl = "response", predict.type)

\method{as.mcmc}{ref.grid}(x, names = TRUE, ...)
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
Type of prediction desired (except in \code{print.xtable}). This only has an effect if there is a known transformation or link function. \code{"response"} specifies that the inverse transformation be applied. Other valid values are \code{"link"}, \code{"lp"}, and \code{"linear"}; these are equivalent, and request that results be shown for the linear predictor. The default is \code{"link"}, unless the \code{"predict.type"} option is in force; see \code{\link{lsm.options}}.

Note that \code{type} is also an argument for the \code{print.xtable} method; it is passed to \code{\link[xtable]{print.xtableList}} in the \pkg{xtable} package.
}
  \item{df}{
If non-missing a constant number of degrees of freedom to use in constructing confidence intervals and \emph{P} values (\code{NA} specifies asymptotic results).
}
  \item{null}{Null hypothesis value(s) against which estimates are tested. May be a single value used for all, or a numeric vector of length equal to the number of tests in each family (i.e., \code{by} group in the displayed table).}
  \item{delta}{Numeric value. If zero, ordinary tests of significance are performed. If positive, this specifies a threshold for testing equivalence (using the TOST or two-one-sided-test method), non-inferiority, or non-superiority, depending on \code{side}. See Details for how the test statistics are defined.}
  \item{side}{Numeric or character value specifying whether the test is left-tailed (\code{-1}, \code{"-"}, code{"<"}, \code{"left"}, or \code{"nonsuperiority"}); right-tailed (\code{1}, \code{"+"}, \code{">"}, \code{"right"}, or \code{"noninferiority"}); or two-sided (\code{0}, \code{2}, \code{"!="}, \code{"two-sided"}, \code{"both"}, \code{"equivalence"}, or \code{"="}).}
  \item{deparse.level}{This argument is needed by the generic \code{rbind} method, but ignored by its \code{ref.grid} method.}
  \item{x}{
The object to be subsetted, printed, plotted, or converted.
}
  \item{y}{This argument is ignored.}
  \item{i}{Integer index(es) of which linear functions to extract.}
  \item{horizontal}{Determines orientation of plotted confidence intervals.}
  \item{intervals}{If \code{TRUE}, confidence intervals are plotted for each estimate}
  \item{comparisons}{If \code{TRUE}, \dQuote{comparison arrows} are added to the plot, in such a way that the degree to which arrows overlap reflects as much as possible the significance of the comparison of the two estimates.}
  \item{alpha}{The \code{alpha} argument to use in constructing comparison arrows.}
  \item{int.adjust}{the multiplicity adjustment method for the plotted confidence intervals; if missing, it defaults to the object's internal \code{adjust} setting (see \code{\link{update}}). (Note: the \code{adjust} argument in \code{plot} sets the adjust method for the comparison arrows, not the confidence intervals.)}
  \item{transform}{Character value. If \code{"response"}, the inverse transformation is applied to the estimates in the grid; if \code{"log"}, the results are formulated as if the response had been \code{log}-transformed; if \code{"none"}, predictions thereof are on the same scale as in \code{object}, and any internal transformation information is preserved. For compatibility with past versions, \code{transform} may also be logical; \code{TRUE} is taken as \code{"response"}, and \code{FALSE} as \code{"none"}.}
  \item{inv.log.lbl}{Character value. This applies only when \code{transform = "log"}, and is used to label the predictions if subsequently summarized with \code{type = "response"}.}
  \item{predict.type}{Character value. If provided, the returned object is first \code{\link{update}}d with the given type, e.g., \code{"response"}.}
  \item{names}{Logical scalar or vector specifying whether variable names are appended to levels in the column labels for the \code{as.mcmc} result -- e.g., column names of \code{treat A} and \code{treat B} versus  just \code{A} and \code{B}. When there is more than one variable involved, the elements of \code{names} are used cyclically.}
\item{\dots, digits, quote, right, caption, label, align, display, auto, include.rownames, sanitize.message.function, xlab, ylab, layout}{For summaries, these are additional arguments passed to other methods including \code{\link{print.data.frame}}, 
\code{\link{xtableList}}, \code{\link{print.xtableList}},
\code{\link{update}}, or \code{\link{dotplot}} as appropriate. If not specified, appropriate defaults are used. For example, the default \code{layout} is one column of horizontal panels or one row of vertical panels.}
}
\details{
\bold{Defaults for summarization, etc.:}
The \code{misc} slot in \code{object} contains default values for \code{by}, \code{infer}, \code{level}, \code{adjust}, \code{type}, \code{null}, \code{side}, and \code{delta}. These defaults vary depending on the code that created the object. The \code{\link{update}} method may be used to change these defaults. In addition, any options set using \samp{lsm.options(summary=...)} will trump those stored in the object's \code{misc} slot.

\bold{Transformations and links:}
With \code{type="response"}, the transformation assumed can be found in \samp{object@misc$tran}, and its label, for the summary is in \samp{object@misc$inv.lbl}. At this time, \code{tran} must be one of the named transformations valid for \code{\link{make.link}}. Any \eqn{t} or \eqn{z} tests are still performed on the scale of the linear predictor, not the inverse-transformed one. Similarly, confidence intervals are computed on the linear-predictor scale, then inverse-transformed.

\bold{Confidence-limit and P-value adjustments:}
The \code{adjust} argument specifies a multiplicity adjustment for tests or confidence intervals. This adjustment always is applied \emph{separately} to each table or sub-table that you see in the printed output (see the details on \code{rbind} below for how to combine tables). The valid values of \code{adjust} are as follows:
\describe{
\item{\code{"tukey"}}{Uses the Studentized range distribution with the number of means in the family. (Available for two-sided cases only.)} 
\item{\code{"scheffe"}}{Computes \eqn{p} values from the \eqn{F} distribution, according to the Scheffe critical value of \eqn{\sqrt{kF(k,d)}}{sqrt[k*F(k,d)]}, where \eqn{d} is the error degrees of freedom and \eqn{k} is (family size minus 1) for contrasts, and (number of estimates) otherwise. (Available for two-sided cases only.)}
\item{\code{"sidak"}}{Makes adjustments as if the estimates were independent (a conservative adjustment in many cases).}
\item{\code{"bonferroni"}}{Multiplies \eqn{p} values, or divides significance levels by the number of estimates. This is a conservative adjustment.}
\item{\code{"dunnettx"}}{Uses an approximation to the Dunnett distribution for a family of estimates having pairwise correlations of \eqn{0.5} (as is true when comparing treatments with a control with equal sample sizes). The accuracy of the approximation improves with the number of simultaneous estimates, and is much faster than \code{"mvt"}. (Available for two-sided cases only.)}
\item{\code{"mvt"}}{Uses the multivariate \eqn{t} distribution to assess the probability or critical value for the maximum of \eqn{k} estimates. This method produces the same \eqn{p} values and intervals as the default \code{summary} or \code{confint} methods to the results of \code{\link{as.glht}}. In the context of pairwise comparisons or comparisons with a control, this produces \dQuote{exact} Tukey or Dunnett adjustments, respectively. However, the algorithm (from the \pkg{mvtnorm} package) uses a Monte Carlo method, so results are not exactly repeatable unless the random-number seed is used (see \code{\link[base]{set.seed}}). As the family size increases, the required computation time will become noticeable or even intolerable, making the \code{"tukey"}, \code{"dunnettx"}, or others more attractive.}
\item{\code{"none"}}{Makes no adjustments to the \eqn{p} values.}
} % end \describe {}

For P-value adjustments only, the Bonferroni-inequality-based adjustment methods in \code{\link{p.adjust}} are also available (currently, these include \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, and \code{"none"}). If a \code{p.adjust.methods} method other than \code{"bonferroni"} or \code{"none"} is specified for confidence limits, the straight Bonferoni adjustment is used instead.
Also, if an adjustment method is not appropriate (e.g., using \code{"tukey"} with one-sided tests, or with results that are not pairwise comparisons), a more appropriate method (usually \code{"sidak"}) is substituted. 

In some cases, confidence and \eqn{p}-value adjustments are only approximate -- especially when the degrees of freedom or standard errors vary greatly within the family of tests. The \code{"mvt"} method is always the correct one-step adjustment, but it can be very slow. One may use \code{\link{as.glht}} with methods in the \pkg{multcomp} package to obtain non-conservative multi-step adjustments to tests.

\bold{Information:}
The \code{str} method outputs a very brief summary of the object, whereas \code{levels} produces a \code{data.frame} containing the combinations of levels of predictor values that define the reference grid.

\bold{\code{xtable}-related methods:}
The \code{\link{xtable}} methods actually use \code{\link{xtableList}}, because of the ability to display messages such as those for P-value adjustments. These methods return an object of class \code{"xtable.lsm"} -- an extension of \code{"xtableList"}. Unlike other \code{xtable} methods, the number of digits defaults to 4; and degrees of freedom and \emph{t} ratios are always formatted independently of \code{digits}. The \code{print} method uses \code{\link{print.xtableList}}, and any \code{\dots} arguments are passed there.

\bold{\code{rbind} and \code{[} methods:}
\code{rbind} can be used to combine two or more reference grids into one. The \code{"["} method for \code{ref.grid}s may be used to obtain a subset. The primary reason for doing this would be to redefine the family of tests to which a P-value adjustment method applies. In \code{rbind}, the variables defined in the objects' grids are merged into one grid, and the returned object has no \dQuote{by} variables and the multiplicity adjustment method set to \code{"mvt"} (as this is likely the only appropriate one).
\code{rbind} throws an error if there are any mismatches among the dimensions, fixed-effect coefficients, or covariance matrices.

\bold{Non-estimable cases:}
When the model is rank-deficient, each row \code{x} of \code{object}'s \code{linfct} slot is each checked for estimability. If \code{sum(x*bhat)} is found to be non-estimable, then an \code{NA} is displayed for the estimate (as well as any associated statistics). This check is performed using the orthonormal basis \code{N} in the \code{nbasis} slot for the null space of the rows of the model matrix. Estimability fails when \eqn{||Nx||^2 / ||x||^2} exceeds \code{tol}, which by default is \code{1e-8}. You may change it via \code{\link{lsm.options}} by setting \code{estble.tol} to the desired value.

\bold{More on tests:}
When \code{delta = 0}, test statistics are of the usual form \samp{(estimate - null)/SE}, or notationally, \eqn{t = (Q - \theta_0)/SE} where \eqn{Q} is our estimate of \eqn{\theta}; then left, right, or two-sided \eqn{p} values are produced. 

When \code{delta} is positive, the test statistic depends on \code{side} as follows. 

Left-sided (nonsuperiority, \eqn{H_0: \theta \ge \theta_0 + \delta} versus \eqn{H_1: \theta < \theta_0 + \delta}): \eqn{t = (Q - \theta_0 - \delta)/SE}. The \eqn{p} value is the lower-tail probability.

Right-sided (noninferiority): \eqn{H_0: \theta \le \theta_0 - \delta} versus \eqn{H_1: \theta > \theta_0 - \delta}): \eqn{t = (Q - \theta_0 + \delta)/SE}. The \eqn{p} value is the upper-tail probability.

Two-sided (equivalence): \eqn{H_0: |\theta - \theta_0| \ge \delta} versus \eqn{H_1: |\theta - \theta_0| < \delta}): \eqn{t = (|Q - \theta_0| - \delta)/SE}. The \eqn{p} value is the \emph{lower}-tail probability.


\bold{Plots:} 
The \code{plot} method for \code{"lsmobj"} or \code{"summary.ref.grid"} objects (but not \code{"ref.grid"} objects themselves) produces a plot displaying confidence intervals for the estimates. If any \code{by} variables are in force, the plot is divided into separate panels. These functions use the \code{\link[lattice]{dotplot}} function, and thus require that the \pkg{lattice} package be installed. For \code{"summary.ref.grid"} objects, the \code{\dots} arguments in \code{plot} are passed \emph{only} to \code{dotplot}, whereas for \code{"lsmobj"} objects, the object is updated using \code{\dots} before summarizing and plotting.

In plots with \code{comparisons = TRUE}, the resulting arrows are only approximate, and in some cases may fail to accurately reflect the pairwise comparisons of the estimates -- especially when estimates having large and small standard errors are intermingled in just the wrong way.

\bold{Re-gridding:}
The \code{regrid} function reparameterizes an existing \code{ref.grid} so that its \code{linfct} slot is the identity matrix and its \code{bhat} slot consists of the estimates at the grid points. If \code{transform} is \code{TRUE}, the inverse transform is applied to the estimates. Outwardly, the \code{summary} after applying \code{regrid} is identical to what it was before (using \samp{type="response"} if \code{transform} is \code{TRUE}). But subsequent contrasts will be conducted on the transformed scale -- which is the reason this function exists. See the example below. In cases where the degrees of freedom depended on the linear function being estimated, the d.f. from the reference grid are saved, and a kind of \dQuote{containment} method is substituted in the returned object whereby the calculated d.f. for a new linear function will be the minimum d.f. among those having nonzero coefficients. This is kind of an \emph{ad hoc} method, and it can over-estimate the degrees of freedom in some cases.

\bold{MCMC samplers:}
When the object's \code{post.beta} slot is non-trivial, \code{as.mcmc} will return an \code{\link[coda]{mcmc}} object that can be summarized or plotted using methods in the \pkg{coda} package. Specifically, \code{post.beta} is transformed by post-multiplying by \code{t(linfct)}, creating a sample from the posterior distribution of LS means.
}

\value{
The \code{summary} method for \code{"ref.grid"} objects returns an object of class \code{"summary.ref.grid"}, which extends \code{"data.frame"}. \code{xtable} returns an object of class \code{"xtable.lsm"}, as explained in details. \code{plot} returns an object of class \code{"trellis"}. \code{vcov} returns the covariance matrix of the product of the object's \code{linfct} and \code{bhat} slots. \code{as.mcmc} returns a \pkg{coda} \code{mcmc} object.
}
%%\references{}
\author{
Russell V. Lenth
}
%\note{}

\seealso{
Methods for the closely related \code{"lsmobj"} class can be found in \code{\link[lsmeans]{contrast}}, \code{\link[lsmeans]{cld}}, and \code{\link[lsmeans]{glht}}. For more on Bonferroni-based P-value adjustments, see \code{\link{p.adjust}}. Also, \code{\link{test}} and \code{\link{confint}} are essentially front-ends for \code{summary}, so additional examples may be found there.
}
\examples{
require(lsmeans)
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.rg <- ref.grid(warp.lm)
str(warp.rg)
levels(warp.rg)

summary(warp.rg)

summary(warp.rg, by = "wool", 
        infer = c(TRUE, FALSE), level = .90, adjust = "sidak")

# Do all pairwise comparisons within rows or within columns, 
# all considered as one faily of tests:
w.t <- pairs(lsmeans(warp.rg, ~ wool | tension))
t.w <- pairs(lsmeans(warp.rg, ~ tension | wool))
rbind(w.t, t.w)

# Transformed response
sqwarp.rg <- ref.grid(update(warp.lm, sqrt(breaks) ~ .))
summary(sqwarp.rg)

# Back-transformed results - compare with summary of 'warp.rg'
summary(sqwarp.rg, type = "response")

# But differences of sqrts can't be back-transformed
summary(pairs(sqwarp.rg, by = "wool"), type = "response")

# We can do it via regrid
sqwarp.rg2 <- regrid(sqwarp.rg)
summary(sqwarp.rg2)  # same as for sqwarp.rg with type = "response"
pairs(sqwarp.rg2, by = "wool")

# Logistic regression
# Reshape the Titanic data
Titan <- do.call("expand.grid", dimnames(Titanic)[-4])
Titan$Died <- matrix(Titanic, ncol=2)
Titan.glm <- glm(Died ~ (Class + Sex + Age)^2, 
    family = binomial, data = Titan)
Titan.lsm <- lsmeans(Titan.glm, ~ Class|Sex, at = list(Age="Adult"))
summary(Titan.lsm, type="response")
summary(pairs(Titan.lsm), type="response")

# Nonsuperiority test: Is any class no more likely to die than
# the 1st class passengers?
summary(contrast(Titan.lsm, "trt.vs.ctrl1"), delta = 1, 
    adjust = "none", side = "<")


# Plot 90% CIs on the response scale
plot(Titan.lsm, type = "response", level = .90, 
     xlab = "Predicted probability of drowning")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }
