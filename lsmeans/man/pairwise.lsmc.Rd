\name{pairwise.lsmc}
\alias{pairwise.lsmc}
\alias{revpairwise.lsmc}
\alias{poly.lsmc}
\alias{trt.vs.ctrl.lsmc}
\alias{trt.vs.ctrl1.lsmc}
\alias{trt.vs.ctrlk.lsmc}

\title{
Contrast families
}
\description{
These functions return standard sets of contrast coefficients.
The name of any of these functions (with the \code{.lsmc} omitted) may be used in the
left-hand side of a \code{spec} formula in \code{\link{lsmeans}}.
}
\usage{
pairwise.lsmc(levs, ...)
revpairwise.lsmc(levs, ...)

poly.lsmc(levs, max.degree = min(6, k - 1))

trt.vs.ctrl.lsmc(levs, ref = 1)
trt.vs.ctrl1.lsmc(levs, ...)
trt.vs.ctrlk.lsmc(levs, ...)
}

\arguments{
  \item{levs}{Vector of factor levels}
  \item{\dots}{Additional arguments, ignored but needed to make these functions interchangeable}
  \item{max.degree}{The maximum degree of the polynomial contrasts in \code{poly.lsmc}}
  \item{ref}{Reference level (or control group) in \code{trt.vs.ctrl.lsmc}}
}

\details{
\code{pairwise.lsmc} and \code{revpairwise.lsmc} generate contrasts for all pairwise comparisons among least-squares means at the levels in \code{levs}. The distinction is in which direction they are subtracted. For factor levels A, B, C, D, \code{pairwise.lsmc} generates the comparisons A-B, A-C, A-D, B-C, B-D, and C-D, wereas \code{revpairwise.lsmc} generates B-A, C-A, C-B, D-A, D-B, and D-C.

\code{poly.lsmc} generates orthogonal polynomial contrasts, assuming equally-spaced factor levels. These are derived from the \code{\link{poly}} function, but an ad hoc algorithm is used to scale them to integer coefficients that are (usually) the same as in published tables of orthogonal polynomial contrasts.

\code{trt.vs.ctrl.lsmc} and its relatives generate contrasts for comparing one level (or the average over specified levels) with each of the other levels. The argument \code{ref} should be the \emph{index}(es) (not the labels) of the reference level(s). \code{trt.vs.ctrl1.lsmc} is the same as \code{trt.vs.ctrl , ref=1)}, and \code{trt.vs.ctrlk.lsmc} is the same as \code{trt.vs.ctrl , ref=length(levs))}
}
\value{
A \code{data.frame}, each column containing contrast coefficients for \code{levs}.
The \code{"desc"} attribute is used to label the results in \code{lsmeans},
and the \code{"adjust"} attribute gives the default adjustment method for multiplicity.
}

\author{
Russell V. Lenth, The University of Iowa
}
\note{
You may create your own contrast functions, using these as guides. A function named \code{mycont.lsmc} may be invoked in \code{lsmeans} via \code{lsmeans(\var{object}, mycontr ~ \var{factor})}. The \code{"desc"} and \code{"adjust"} attributes are optional; if absent, the root name of the function is used as \code{"desc"}, and no adjustment is requested for p values. See the examples.
}


\seealso{
\code{\link{lsmeans}}
}
\examples{
### View orthogonal polynomials for 4 levels
poly.lsmc(1:4)

\dontrun{
### Setting up a custom contrast function
helmert.lsmc <- function(levs, ...) {
  M <- as.data.frame(contr.helmert(levs))
  names(M) <- paste(levs[-1],"vs earlier")
  attr(M, "desc") <- "Helmert contrasts"
  M
}
lsmeans(Oats.lme, helmert ~ Variety)
}

}

\keyword{ models }
\keyword{ regression }
\keyword{ htest }
