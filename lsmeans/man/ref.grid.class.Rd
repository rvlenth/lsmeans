\name{ref.grid-class}
\Rdversion{1.1}
\docType{class}
\alias{ref.grid-class}
\alias{lsmobj-class}
\alias{lsmeans,ref.grid,character-method}
\alias{show,ref.grid-method}
\alias{show,lsmobj-method}
\alias{summary,ref.grid-method}

\title{Classes \code{"ref.grid"} and \code{"lsmobj"}}
\description{
A reference grid encapsulates everything needed to compute least-squares means, independently of the underlying model object. The \code{"lsmobj"} class is a minor extension of \code{"ref.grid"} where the linear predictors for the reference grid are transformed in some linear way such as marginal averages or contrasts.
}

\section{Objects from the Classes}{
Objects of class \code{"ref.grid"} are most commonly created by calling the \code{\link{ref.grid}} function. 

Objects of class \code{"lsmobj"} are created by calling \code{\link{lsmeans}} or a related function such as \code{\link{contrast}}.}

\section{Slots}{
  \describe{
    \item{\code{model.info}:}{Object of class \code{"list"} containing the elements \code{call} (the call that produced the model), \code{terms} (its \code{terms} object), and \code{xlev} (factor-level information)}
    \item{\code{roles}:}{Object of class \code{"list"} containing at least the elements \code{predictors}, \code{responses}, and \code{multresp}. These are character vectors of names of these variables.}
    \item{\code{grid}:}{Object of class \code{"data.frame"} containing the combinations of the variables that define the reference grid. In addition, there is an auxiliary column named \code{".freq."} holding the observed frequencies of each factor combination (excluding covariates). If the model has one or more \code{\link{offset}()} calls, there is an another auxiliary column named \code{".offset."}. Auxiliary columns are not considered part of the reference grid. (However, any variables included in \code{offset} calls \emph{are} in the reference grid.)}
    \item{\code{levels}:}{Object of class \code{"list"} with each entry containing the distinct levels of variables in the reference grid. Note that \code{grid} is obtained by applying the function \code{\link{expand.grid}} to this list}
    \item{\code{matlevs}:}{Object of class \code{"list"} Like \code{levels} but has the levels of any matrices in the original dataset. Matrix columns must always be reduced to a single value for purposes of the reference grid }
    \item{\code{linfct}:}{Object of class \code{"matrix"} giving the linear functions of the regression coefficients for predicting each element of the reference grid. The rows of this matrix go in one-to-one correspondence with the rows of \code{grid}, and the columns with elements of \code{bhat}}
    \item{\code{bhat}:}{Object of class \code{"numeric"} with the regression coefficients. If there is a multivariate response, this must be flattened to a single vector, and \code{linfct} and \code{V} redefined appropriately. Important: \code{bhat} must \emph{include} any \code{NA} values produced by collinearity in the predictors. These are taken care of later in the estimability check.}
    \item{\code{nbasis}:}{Object of class \code{"matrix"} with the basis for the non-estimable functions of the regression coefficients. Every LS mean will correspond to a linear combination of rows of \code{linfct}, and that result must be orthogonal to all the columns of \code{nbasis} in order to be estimable. This will be \code{NULL} if everything is estimable}
    \item{\code{V}:}{Object of class \code{"matrix"}, the symmetric variance-covariance matrix of \code{bhat} }
    \item{\code{dffun, dfargs:}}{Objects of class \code{"function"} and \code{"list"} respectively. \code{dffun(k,dfargs)} should return the degrees of freedom for the linear function \code{sum(k*bhat)}, or \code{NA} if unavailable}
    \item{\code{misc}:}{A \code{list} containing additional information used by methods. These include at least the following: \code{estName} (the label for the estimates of linear functions), and the default values of \code{infer}, \code{level}, and \code{adjust} to be used in the \code{\link{summary}} method. Elements in this slot may be modified if desired using the \code{\link{update}} method.}
}} % end of describe and slots section

\section{Extends}{
\code{Class "lsmobj"} extends \code{Class "ref.grid"}, directly. There is hardly a difference between these classes except for how the slots \code{linfct} and \code{grid} are obtained, and their \code{show} methods.}

\section{Methods}{
  All methods for these objects are S3 methods except for \code{show}.
  \describe{
    \item{\code{show}:}{Prints the results of \code{str} for \code{ref.grid} objects, and \code{summary} for \code{lsmobj} objects.}
    \item{\code{str}:}{Displays a brief listing of the variables and levels defining the grid. }
    \item{\code{summary}:}{Displays a summary of estimates, standard errors, degrees of freedom, and optionally, tests and/or confidence intervals. }
    \item{\code{lsmeans}:}{Computes least-squares means and creates an \code{"lsmobj"} object.}
    \item{\code{confint}:}{Confidence intervals for lsmeans.} 
    \item{\code{test}:}{Hypothesis tests. }
    \item{\code{cld}:}{Compact-letter display for tests of pairwise comparisons}
    \item{\code{contrast}:}{Contrasts among lsmeans. }
    \item{\code{pairs}:}{A special case of \code{contrasts} for pairwise comparisons. }
    \item{\code{update}:}{Change defaults used primarily by \code{summary}, such as transformation, p-value adjustment, and confidence level.}
} % end of \describe
} % end of Methods section


%\references{}
\author{
Russell V. Lenth
}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{ref.grid}}, \code{\link{lsmeans}}
}

\examples{
showClass("ref.grid")
showClass("lsmobj")
}
\keyword{classes}
