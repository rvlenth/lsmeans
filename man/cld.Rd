\name{cld}
\alias{cld}
\alias{cld.ref.grid}

\title{Compact letter display of pairwise comparisons}

\description{
Extract and display information on all pairwise comparisons of least-squares means.}

\usage{
\method{cld}{ref.grid}(object, details = FALSE, sort = TRUE, by, alpha = 0.05, 
    Letters = c("1234567890", LETTERS, letters), ...)
}

\arguments{
  \item{object}{An object of class \code{ref.grid}}
  \item{details}{Logical value determining whether detailed
  information on tests of pairwise comparisons is displayed}
  \item{sort}{Logical value determining whether the LS means are 
  sorted before the comparisons are produced}
  \item{by}{Character value giving the name or names of variables by which separate  
  families of comparisons are tested. 
  If \code{NULL}, all means are compared. If missing, and a \code{by} variable was used in creating \code{object}, it is used as the \code{by} variable in \code{cld}.}
  \item{alpha}{Numeric value giving the significance level for the comparisons}
  \item{Letters}{Character vector of letters to use in the display.
  Any strings of length greater than 1 are expanded into individual characters}
  \item{\dots}{Arguments passed to \code{\link{contrast}} (for example, an \code{adjust} method)}
}
\details{
This function uses the Piepho (2004) algorithm (as implemented in the \pkg{multcompView} package) to generate a compact letter display of all pairwise comparisons of least-squares means. The function obtains (possibly adjusted) \emph{P} values for all pairwise comparisons of means, using the \code{\link{contrast}} function with \code{method = "pairwise"}. When a \code{P} value exceeds \code{alpha}, then the two means have at least one letter in common.
}
\value{
When \code{details == FALSE}, an object of class \code{summary.ref.grid} (which inherits from \code{data.frame}) showing the summary of LS means with an added column named \code{.groups} with the cld information. When \code{details == TRUE}, a list the object just described, as well as the summary of the \code{contrast} results showing each comparison, its estimate, standard error, \eqn{t} ratio, and adjusted \eqn{P} value.
}
\references{
Hans-Peter Piepho (2004) An algorithm for a letter-based representation of all pairwise comparisons, \emph{Journal of Computational and Graphical Statistics}, 13(2), 456-466.
}
\author{Russell V. Lenth}
\note{This function requires the \pkg{multcompView} package to be installed. Otherwise an error message is produced.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[multcomp]{cld}} in the \pkg{multcomp} package
}
\examples{
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.lsm <- lsmeans(warp.lm, ~ tension | wool)
cld(warp.lsm)                  # implicitly uses by = "wool"
cld(warp.lsm, by = "tension")  # overrides implicit 'by'

# Mimic grouping bars and compare all 6 means
cld(warp.lsm, by = NULL, Letters = "||||||||", alpha = .01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
