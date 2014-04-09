\name{update.ref.grid}
\alias{update.ref.grid}
\alias{update}

\title{
Update a \code{ref.grid} or \code{lsmobj} object
}
\description{
Objects of class \code{ref.grid} or \code{lsmobj} contain several settings in their \code{"misc"} slot that affect primarily the 
defaults used by \code{\link{summary}}. This \code{update} method allows them to be changed more safely than by modifying this slot directly.
}
\usage{
\S3method{update}{ref.grid}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{ref.grid} (or its extension, \code{lsmobj})
}
  \item{\dots}{
Arguments specifying elements' names and their new values.
}
}
\details{
The names in \code{\dots} are partially matched against those that are valid, and if a match is found, it adds or replaces the current setting. The valid names are

\describe{
\item{tran}{(\code{list} or \code{character}) specifies the transformation which, when inverted, determines the results displayed by \code{\link{summary}}, \code{\link{predict}}, or \code{\link{lsmip}} when \code{type="response"}. The value may be the name of a standard transformation from \code{\link{make.link}}, or, for a custom transformation, a \code{list} containing at least the functions \code{linkinv} (the inverse of the transformation) and \code{mu.eta} (the derivative thereof). See the Examples. }

\item{estName}{(\code{character}) is the column label used for displaying predictions or LS means.}

\item{inv.label}{(\code{character)}) is the column label to use for predictions or LS means when \code{type="response"}.}

\item{by.vars}{(\code{character)} vector or \code{NULL}) the variables used for grouping in the summary, and also for defining subfamilies in a call to \code{\link{contrast}}.}

\item{pri.vars}{(\code{character} vector) are the names of the grid variables that are not in \code{by.vars}. Thus, the combinations of their levels are used as columns in each table produced by \code{\link{summary}}.}

\item{adjust}{(\code{character)}) is the default for the \code{adjust} argument in \code{\link{summary}}.}

\item{famSize}{(integer) is the \code{nmeans} parameter for \code{\link{ptukey}} when \code{adjust="tukey"}. }

\item{infer}{(\code{logical} vector of length 2) is the default value of \code{infer} in \code{\link{summary}}.}

\item{level}{(numeric) is the default confidence level, \code{level}, in \code{\link{summary}}}

\item{avgd.over}{(\code{character)} vector) are the names of the variables whose levels are averaged over in obtaining marginal averages of predictions, i.e., LS means. Changing this might produce a misleading printout, but setting it to \code{character(0)} will suppress the \dQuote{averaged over} message in the summary.}

\item{methDesc}{(\code{character}) is a string that may used for creating names for a list of \code{lsmobj} objects. }

} % end \describe
} % end details
\value{
A copy of \code{object} with its \code{"misc"} slot modified.
}

\author{
Russell V. Lenth
}
\note{Unlike the \code{update} method for model classes (\code{lm}, \code{glm}, etc.), this does not re-fit or re-estimate anything; but it does affect how \code{object} is treated by other methods for its class.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary}}
}
\examples{
# An altered log transformation

warp.lm1 <- lm(log(breaks + 1) ~ wool*tension, data = warpbreaks)
rg1 <- update(ref.grid(warp.lm1), 
        tran = list(linkinv = function(eta) exp(eta) - 1,
                    mu.eta = function(eta) exp(eta)),
        inv.lbl = "pred.breaks")
        
summary(rg1, type = "response")        
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }% __ONLY ONE__ keyword per line
