\name{pairwise.lsmc}
\alias{pairwise.lsmc}
\alias{revpairwise.lsmc}
\alias{poly.lsmc}
\alias{trt.vs.ctrl.lsmc}
\alias{trt.vs.ctrl1.lsmc}
\alias{trt.vs.ctrlk.lsmc}
\alias{eff.lsmc}
\alias{del.eff.lsmc}
\alias{tukey.lsmc}
\alias{dunnett.lsmc}
\alias{consec.lsmc}
\alias{mean_chg.lsmc}

\title{
Contrast families
}
\description{
These functions return standard sets of contrast coefficients.
The name of any of these functions (with the \code{.lsmc} omitted) may be used as the \code{method} argument in \code{\link{contrast}}, or as the \code{contr} argument or left-hand side of a \code{spec} formula in \code{\link{lsmeans}}.
}
\usage{
pairwise.lsmc(levs, ...)
revpairwise.lsmc(levs, ...)
tukey.lsmc(levs, reverse = FALSE)

poly.lsmc(levs, max.degree = min(6, k - 1))

trt.vs.ctrl.lsmc(levs, ref = 1)
trt.vs.ctrl1.lsmc(levs, ...)
trt.vs.ctrlk.lsmc(levs, ...)
dunnett.lsmc(levs, ref = 1)

consec.lsmc(levs, reverse = FALSE, ...)
mean_chg.lsmc(levs, reverse = FALSE, ...)

eff.lsmc(levs, ...)
del.eff.lsmc(levs, ...)
}

\arguments{
  \item{levs}{Vector of factor levels}
  \item{\dots}{Additional arguments, ignored but needed to make these functions interchangeable}
  \item{max.degree}{The maximum degree of the polynomial contrasts in \code{poly.lsmc}}
  \item{reverse}{Logical value to determine the direction of comparisons, e.g., pairwise (if \code{TRUE}) or reverse-pairwise (if \code{FALSE}) comparisons.}
  \item{ref}{Reference level (or control group) in \code{trt.vs.ctrl.lsmc}}
}

\details{
Each contrast family has a default multiple-testing adjustment as noted below. These adjustments are often only approximate; for a more exacting adjustment, use the interfaces provided to \code{\link[multcomp]{glht}} in the \pkg{multcomp} package.

\code{pairwise.lsmc}, \code{revpairwise.lsmc}, and \code{tukey.lsmc} generate contrasts for all pairwise comparisons among least-squares means at the levels in \code{levs}. The distinction is in which direction they are subtracted. For factor levels A, B, C, D, \code{pairwise.lsmc} generates the comparisons A-B, A-C, A-D, B-C, B-D, and C-D, whereas \code{revpairwise.lsmc} generates B-A, C-A, C-B, D-A, D-B, and D-C. \code{tukey.lsmc} invokes \code{pairwise.lsmc} or \code{revpairwise.lsmc} depending on \code{reverse}. The default multiplicity adjustment method is \code{"tukey"}, which is approximate when the standard errors differ.

\code{poly.lsmc} generates orthogonal polynomial contrasts, assuming equally-spaced factor levels. These are derived from the \code{\link{poly}} function, but an ad hoc algorithm is used to scale them to integer coefficients that are (usually) the same as in published tables of orthogonal polynomial contrasts. The default multiplicity adjustment method is \code{"none"}.

\code{trt.vs.ctrl.lsmc} and its relatives generate contrasts for comparing one level (or the average over specified levels) with each of the other levels. The argument \code{ref} should be the \emph{index}(es) (not the labels) of the reference level(s). \code{trt.vs.ctrl1.lsmc} is the same as \code{trt.vs.ctrl} with a reference value of \code{1}, and \code{trt.vs.ctrlk.lsmc} is the same as \code{trt.vs.ctrl} with a reference value of \code{length(levs)}. \code{dunnett.lsmc} is the same as \code{trt.vs.ctrl}.
The default multiplicity adjustment method is \code{"dunnettx"}, a close approximation to the Dunnett adjustment.

\code{consec.lsmc} and \code{mean_chg.lsmc} are useful for contrasting treatments that occur in sequence. For a factor with levels A, B, C, D, E, \code{consec.lsmc} generates the comparisons B-A, C-B, and D-C, while \code{mean_chg.lsmc} generates the contrasts (B+C+D)/3 - A, (C+D)/2 - (A+B)/2, and D - (A+B+C)/3.  With \code{reverse = TRUE}, these differences go in the opposite direction.

\code{eff.lsmc} and \code{del.eff.lsmc} generate contrasts that compare each level with the average over all levels (in \code{eff.lsmc}) or over all other levels (in \code{del.eff.lsmc}). These differ only in how they are scaled. For a set of \eqn{k} lsmeans, \code{del.eff.lsmc} gives weight \eqn{1} to one lsmean and weight \eqn{-1/(k-1)} to the others, while \code{eff.lsmc} gives weights \eqn{(k-1)/k} and \eqn{-1/k} respectively, as in subtracting the overall lsmean from each lsmean.
The default multiplicity adjustment method is \code{"fdr"}. This is a Bonferroni-based method and is slightly conservative; see \code{\link{p.adjust}}
}
\value{
A \code{data.frame}, each column containing contrast coefficients for \code{levs}.
The \code{"desc"} attribute is used to label the results in \code{lsmeans},
and the \code{"adjust"} attribute gives the default adjustment method for multiplicity.
}

\author{
Russell V. Lenth
}
\note{
You may create your own contrast functions, using these as guides. A function named \code{mycontr.lsmc} may be invoked in \code{lsmeans} via, e.g., \preformatted{lsmeans(\var{object}, mycontr ~ \var{factor})} 
The \code{"desc"}, \code{"adjust"}, and \code{"offset"} attributes are optional; if present, these are passed to \code{contrast}. If absent, the root name of the function is used as \code{"desc"}, and no adjustment is requested for p values. See the examples.
}


\seealso{
\code{\link{lsmeans}}, \code{\link[multcomp]{glht}}
}
\examples{
### View orthogonal polynomials for 4 levels
poly.lsmc(1:4)

### Setting up a custom contrast function
helmert.lsmc <- function(levs, ...) {
  M <- as.data.frame(contr.helmert(levs))
  names(M) <- paste(levs[-1],"vs earlier")
  attr(M, "desc") <- "Helmert contrasts"
  M
}
warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
lsmeans(warp.lm, helmert ~ tension | wool)
}

\keyword{ models }
\keyword{ regression }
\keyword{ htest }
