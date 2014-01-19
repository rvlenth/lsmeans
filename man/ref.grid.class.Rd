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
A reference grid encapsulates everything needed to compute least-squares means, independently of the underlying model object. The \code{"lsmobj"} class is a minor extension of \code{"ref.grid"} where the linear predictors for the reference grid are averaged over the levels of zero or more variables in the grid.
}

\section{Objects from the Classes}{
Objects of class \code{"ref.grid"} are most commonly created by calling the \code{\link{ref.grid}} function. 

Objects of class \code{"lsmobj"} are created by calling \code{\ref{lsmeans}} or a related function such as \code{\ref{contrasts}}.}

\section{Slots}{
  \describe{
    \item{\code{model.info}:}{Object of class \code{"list"} containing the elements \code{call} (the call that produced the model), \code{terms} (its \code{terms} object), and \code{xlev} (factor-level information)}
    \item{\code{roles}:}{Object of class \code{"list"} containing at least the elements \code{predictors}, \code{responses}, and \code{multresp}. These are character vectors of names of these variables.}
    \item{\code{grid}:}{Object of class \code{"data.frame"} containing the combinations of the variables that define the reference grid.}
    \item{\code{levels}:}{Object of class \code{"list"} with each entry containing the distinct levels of variables in the reference grid. Note that \code{grid} is obtained by applying the function \code{\ref{expand.grid}} to this list}
    \item{\code{matlevs}:}{Object of class \code{"list"} Like \code{levels} but has the levels of any matrices in the original dataset. Matrix columns must always be reduced to a single value for purposes of the reference grid }
    \item{\code{linfct}:}{Object of class \code{"matrix"} giving the linear functions of the regression coefficients for predicting each element of the reference grid. The rows of this matrix go in one-to-one correspondence with the rows of \code{grid}, and the columns with elements of \code{bhat}}
    \item{\code{bhat}:}{Object of class \code{"numeric"} with the regression coefficients. If there is a multivariate response, this must be flattened to a single vector, and \code{linfct} and \code{V} redefined appropriately. Important: \code{bhat} must \emph{include} any \code{NA} values produced by collinearity in the predictors. These are taken care of later in the estimability check.}
    \item{\code{nbasis}:}{Object of class \code{"matrix"} with the basis for the non-estimable functions of the regression coefficients. Every LS mean will correspond to a linear combination of rows of \code{linfct}, and it must be orthogonal to all the columns of \code{nbasis} in order to be estimable. This will be \code{NULL} if everything is estimable}
    \item{\code{V}:}{Object of class \code{"matrix"}, the symmetric variance-covariance matrix of \code{bhat} }
    \item{\code{ddfm}:}{Object of class \code{"function"}. \code{ddfm(k,misc)} should return the degrees of freedom for the linear function \code{sum(k*bhat)}, or \code{NA} if unavailable}
    \item{\code{misc}:}{Object of class \code{"list"} passed to \code{ddfm} in case additional arguments are needed. It also serves as a repository of additional information used by methods. These include at least the following: \code{estName} (the label for the estimates of linear functions), and the values of \code{infer}, \code{level}, and \code{adjust} used in the \code{summary} method.}

\section{Extends}{
Class \code{"lsmobj"}} extends \code{Class \code{"ref.grid"}}, directly. There is hardly a difference between these classes except how you think about them.}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ref.grid")}: Displays a brief listing of the variables and levels defining the grid. }
    \item{summary}{\code{signature(object = "ref.grid")}: Displays a summary of estimates, standard errors, degrees of freedom, and optionally, tests and/or confidence intervals. }
    \item{lsmeans}{\code{signature(object = "ref.grid", specs = "character")}: Computes least-squares means and creates an \code{"lsmobj"} object.}
    \item{show}{\code{signature(object = "lsmobj")}: This is the same as \code{summary} for an \code{lsmobj} object. }
    
    %%% THERE WILL BE MORE -- contrasts, tests - convert to S4 methods
    %%% Also mention S3 methods confint, pairs

} % end of \describe
} % end of Methods section


%\references{}
\author{
Russell V. Lenth
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ref.grid}}, \code{\link{lsmeans}}
}

\examples{
showClass("ref.grid")
}
\keyword{classes}
