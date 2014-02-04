\name{MOats}
\alias{MOats}
\docType{data}
\title{
Oats data in multivariate form
}
\description{
This is the \code{Oats} dataset provided in the \pkg{nlme} package,
but it is rearranged as one multivariate observation per plot.
}
\usage{data(MOats)}
\format{
  A data frame with 18 observations on the following 3 variables.
  \describe{
    \item{\code{Variety}}{a factor with levels \code{Golden Rain}, \code{Marvellous}, \code{Victory}}
    \item{\code{Block}}{an ordered factor with levels \code{VI} < \code{V} < \code{III} < \code{IV} < \code{II} < \code{I}}
    \item{\code{yield}}{a matrix with 4 columns, giving the yields with nitrogen concentrations of 0, .2, .4, and .6.}
  }
}
\details{
These data arise from a split-plot experiment reported by Yates (1935) and used as an example in Pinheiro and Bates (2000) and other texts. Six blocks were divided into three whole plots, randomly assigned to the three varieties of oats. The whole plots were each divided into 4 split plots and randomized to the four concentrations of nitrogen.
}
\source{
The dataset \code{\link[nlme]{Oats}} in the \pkg{nlme} package.
}
\references{
Pinheiro, J. C. and Bates D. M. (2000) \emph{Mixed-Effects Models in S and S-PLUS}, Springer, New York. (Appendix A.15)

Yates, F. (1935) Complex experiments, \emph{Journal of the Royal Statistical Society} Suppl. 2, 181-247
}
\examples{
require(lsmeans)
MOats.lm <- lm (yield ~ Block + Variety, data = MOats)
MOats.rg <- ref.grid (MOats.lm, 
                mult.levs = list(nitro = c(0,.2,.4,.6)))
lsmeans(MOats.rg, ~ nitro | Variety)
}
\keyword{datasets}
