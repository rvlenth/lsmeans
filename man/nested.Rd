%    Copyright (c) 2012-2016 Russell V. Lenth                                %

\name{nested}
\alias{nested}

\title{
Efficient factor coding for nested structures
}
\description{

}
\usage{
nested(factor, ..., SAS = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{factor}{The factor that is nested in \code{\dots}.}
  \item{\dots}{The other factor(s) in which \code{factor} is nested.}
  \item{SAS}{Logical value that determines whether the first (when \code{FALSE}) 
      or the last (when \code{TRUE}) level in each nest is regarded as the reference level.}
}
\details{
Any contrast coding (see \code{\link{contrasts}}) will work correctly with either choice of \code{SAS}, in that equivalent models with the same fitted values are produced.
However, for sake of interpretability of regression coefficients, choosing \code{SAS = FALSE}
is compatible with \code{\link{contr.treatment}} contrast coding, 
while \code{SAS = TRUE} is compatible with \code{\link{contr.SAS}} coding. 
} % end of \details

\value{
After determining which levels of \code{factor} occur within each nest, the function creates a factor that has a \eqn{n_i - 1} levels for the \eqn{i}th nest, and a common reference level that is labeled \code{"ref"}. The reference level is lexically the first(last) level of the \code{factor} levels within each nest, depending on whether \code{SAS} is \code{TRUE}(\code{FALSE}).
}

%%%\references{}

\author{
Russell V. Lenth
}

\note{
The use of \code{nested} is needed only when fitting models with nested factors as fixed effects. It is \emph{not} necessary to use \code{nested} in fitting a mixed model, where the nested structure is modeled as a random effect -- via, e.g., \code{lmer} in the \pkg{lme4} package. 

Currently, it is unsafe to include a call to \code{nested} \emph{within} a model formula.
The user is advised to use \code{nested} to create factors for each nested term to be included in a 
model; then fit the model using the names of those created factors. See the example below.
}

\seealso{
  \code{\link{contrasts}}
}
\examples{
require(lsmeans)

# A very small two-period crossover design with 3 subjects one sequence, 2 in the other:
foo <- data.frame(
    seq = factor(c(rep("AB", 6), rep("BA", 4))),
    subj = factor(rep(1:5, each = 2)),
    per = factor(rep(1:2, 5)),
    trt = factor(c(rep(c("A","B"), 3), rep(c("B","A"), 2)))
)
foo <- transform (foo,
    y = c(17.9, 16.2, 17, 15.4, 19.7, 17.9, 19, 20.6, 18.2, 20.4),
    subj.in.seq = nested(subj, seq) )

# Note that 'subj.in.seq' has a ref level plus levels for subjects 2, 3, and 5
foo$subj.in.seq

# Fit a model using R's factory provisions for nested models
foo.lm1 <- lm(y ~ seq/subj + per + trt, data = foo)
coef(foo.lm1)            # rank-deficient model
lsmeans(foo.lm1, pairwise ~ trt)  # LS means not estimable, but difference is

# Fit another model using 'subj.in.seq' instead
foo.lm2 <- lm(y ~ seq + subj.in.seq + per + trt, data = foo)
coef(foo.lm2)            # NOT rank-deficient
lsmeans(foo.lm2, pairwise ~ trt)  # Estimable LS means, same result for difference


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
