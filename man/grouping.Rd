\name{add_grouping}
\alias{add_grouping}


\title{
Create a grouping factor
}
\description{
A new factor \code{G} can be defined having levels based on those of an existing factor \code{A}. This adds \code{G} to the reference grid and redefines \code{A} as being nested in \code{G}.
}
\usage{
add_grouping(object, newname, refname, newlevs)
}
\arguments{
  \item{object}{A \code{\link{ref.grid}} object}
  \item{newname}{Character name of the grouping factor to be created}
  \item{refname}{Character name of an existing factor in \code{object} whose levels are to be grouped}
  \item{newlevs}{Character vector or factor of the same length as the levels for \code{refname}. Each element specifies the corresponding level for \code{newname}}
}
\details{
This function is useful when one wants to group certain levels of some treatment factor together and summarize by those groups -- without fitting a new model. For example, suppose a factor \code{hosp} refers to five hospitals. A call like 
\code{newgrid <- add_grouping(refgrid, "envir", "hosp", c("rural", "urban", "urban", "rural", "urban"))} would add a two-level grouping factor named \code{envir} in which the first and fourth hospitals are rural and the others are urban. Subsequently, we may use \code{lsmeans(newgrid, "envir")} to obtain the marginal means for the rural and urban environments.

If \code{newlevs} is specified as a factor, the user can use its \code{levels} attribute to control the order in which levels are presented. Otherwise, it will be alphabetical.
}
\value{A \code{ref.grid} object with the factor \code{newname} added. It will include a nesting structure such that \code{refname} is nested in \code{newname}, and with \code{newname} added to any existing nesting that involves \code{refname}.}

\author{
Russell V. Lenth
}
\note{If \code{object} already has \code{refname} nested in some other factors, this creates a conflict and therefore it is not allowed; an error is returned.
}



\seealso{\code{\link{ref.grid}}}

\examples{
require("lsmeans")

fiber.lm <- lm(strength ~ diameter + machine, data = fiber)
frg <- ref.grid(fiber.lm)

## Suppose the machines are of two types
gfrg <- add_grouping(frg, newname = "type", refname = "machine", 
                     newlevs = c("computerized", "computerized", "manual"))
gfrg

lsmeans(gfrg, "type", contr = "pairwise")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
