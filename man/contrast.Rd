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
\method{contrast}{ref.grid}(object, method = "eff", interaction = FALSE, 
    by, adjust, offset = NULL, name = "contrast", 
    options = getOption("lsmeans")$contrast, ...)
\method{contrast}{lsm.list}(object, ..., which = 1)

\method{test}{ref.grid}(object, null = 0, joint = FALSE,
    verbose = FALSE, rows, by,  ...)

\method{confint}{ref.grid}(object, parm, level = 0.95, ...)

\method{pairs}{ref.grid}(x, reverse = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{
An object of class \code{"ref.grid"} or its extension, \code{"lsmobj"}.
}
  \item{method}{
Character value giving the root name of a contast method (e.g. \code{"pairwise"}). Alternatively, a named list of contrast coefficients that must each conform to the number of least-squares means in each \code{by} group. This is just like the \code{contr} argument in \code{\link{lsmeans}}. To identify the available methods, see \preformatted{ls("package:lsmeans", pat=".lsmc")} You may define your own \code{.lsmc} function and use its root name as \code{method}. If \code{interaction} is of character type, this argument is ignored.
}
  \item{interaction}{
Character vector or logical value. In multi-factor situations with \code{interaction = FALSE}, the factor combinations are treated as levels of a single \dQuote{uber-factor}, and the contrast specified in \code{method} is applied to it. Otherwise, interaction contrasts are computed: Contrasts are generated for each factor separately, one at a time; and these contrasts are applied to the object (the first time around) or to the previous result (subsequently). (Any factors specified in \code{by} are skipped.) The final result comprises contrasts of contrasts, or, equivalently, products of contrasts for the factors involved. Processing is done in the order of appearance in \code{object@levels}. With \code{interaction = TRUE}, \code{method} (if specified as character) is used for each contrast. If \code{interaction} is a character vector, the elements specify the respective contrast method(s); they are recycled as needed.
}
  \item{by}{
Character names of variable(s) to be used for ``by'' groups. The contrasts or joint tests will be evaluated separately for each combination of these variables. If \code{object} was created with by groups, those are used unless overridden. Use \code{by = NULL} to use no by groups at all.
}
  \item{adjust}{
Method to use for adjusting \emph{P} values. This is passed to \code{\link[lsmeans]{summary}}.
}
  \item{offset}{Numeric vector of the same length as each \code{by} group. These values are added to their respective linear estimates. This argument is ignored when\code{interaction} is not \code{FALSE}.
}
  \item{name}{Name to use to label the contrasts in table headings
  or subsequent contrasts of the returned object. This argument is ignored 
  when\code{interaction} is not \code{FALSE}.
}
  \item{options}{If non-\code{NULL}, a named \code{list} of arguments to pass to \code{\link{update}}, just after the object is constructed.}

  \item{joint}{Logical value. If \code{FALSE}, the arguments are passed to \code{\link{summary}} with \code{infer=c(FALSE,TRUE)}. If \code{TRUE}, a joint test of the hypothesis L beta = null is performed, where L is \code{object@linfct} and beta is the vector of fixed effects estimated by \code{object@betahat}. This will be either an \emph{F} test or a chi-square (Wald) test depending on whether degrees of freedom are available.}
  \item{rows}{Integer values. The rows of L to be tested in the joint test. If missing, all rows of L are used. If not missing, \code{by} variables are ignored.}
  \item{null}{Numeric value specifying the null value(s) being tested against. It may be either a single value, in which case it is used as the null value for all linear functions under test; or a numeric vector of length equal to the number of linear functions.}
  \item{parm}{This is ignored, but it is a required argument of the generic  \code{confint} method.)}
  \item{verbose}{Logical value. If \code{TRUE} and \code{joint==TRUE}, a table of the effects being tested is printed.}
  \item{level}{
Numeric value of the desired confidence level.}
  \item{which}{When \code{object} is a list of \code{lsmobj} objects, this specifies which member of the list is analyzed.}
  \item{reverse}{Logical value determining whether \code{"pairwise"} or \code{"revpairwise"} pairwise comparisons are generated.}
  \item{\dots}{
Additional arguments passed to \code{\link[lsmeans]{summary}} or to a contrast function.}
}
\details{
Though \code{contrast} is ordinarily used to create true contrasts (whose coefficients sum to zero), it may be used to estimate any linear function of the LS means; and \code{offset} expands this capability further by allowing additive constants.
\code{pairs} is equivalent to \code{contrast} with \code{method = "pairwise"}.

\code{confint} and \code{test} (when \code{JOINT==FALSE}) are equivalent to calling \code{\link[lsmeans]{summary}} with \code{infer=c(TRUE,FALSE)} and \code{infer=c(FALSE,TRUE)}, respectively.

When using \code{test} to do a joint test of L beta = null, an error is thrown if any row of L is non-estimable. It is permissible for the rows of L to be linearly dependent as long as \code{null == 0}; a reduced set of contrasts is tested. Linear dependence and nonzero \code{null} cause an error.
}
\value{
\code{contrast} and \code{pairs} return an object of class \code{"lsmobj"}, which is an extension of \code{"ref.grid"}. Consequently, they may be used as arguments to other \code{"lsmobj"} or \code{"ref.grid"} methods. The user may, for example, compute contrasts of contrasts, or re-summarize a set of confidence intervals with a different \code{by} grouping or confidence level.
The ``grid'' for the returned value is simply the set of variables that identify the results. For example, \code{contrast}'s return value is a reference grid for one factor named \code{contrast}.

\code{confint} and \code{test} (when \code{Joint==FALSE}) return an object of class \code{summary.ref.grid}. When \code{JOINT==TRUE}, \code{test} returns a \code{"summary.ref.grid"} object (extends \code{"data.frame"})  with the test statistic, degrees of freedom, and \emph{P} value for each \code{by} group.
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
(warp.pl <- contrast(warp.lsm, "poly", name = "order"))
# Same results with a different adjustment
summary(warp.pl, adjust = "fdr")

# Jointly test the tension effects for each wool
test(warp.pl, joint = TRUE, by = "wool")

# Compare the two contrasts for each order
contrast(warp.pl, "revpairwise", by = "order")

# User-provided contrasts, ignoring the previous by grouping
contrast(warp.lsm, 
    list(c1=c(1,0,0,-1,0,0), c2=c(1,1,1,-1,-1,-1)/3), 
    by = NULL)

# Compare consecutive tension*wool comb's as treatment with 6 levels
contrast(warp.lsm, "consec", by = NULL)

# Interaction contrasts (comparisons of linear and quadratic contrasts)
contrast(warp.lsm, interaction = c("poly", "consec"), by = NULL)
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
\keyword{ htest }

