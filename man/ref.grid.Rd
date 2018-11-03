%    Copyright (c) 2012-2018 Russell V. Lenth                                %

\name{ref.grid}
\alias{ref.grid}
\alias{recover.data}
\alias{lsm.basis}

% others needed to keep checking happy
\alias{summary.ref.grid}
\alias{lsmeans}
\alias{contrast}


\title{
Create a reference grid from a fitted model
}
\description{
These functions are provided in \pkg{lsmeans} because they have been renamed in \pkg{emmeans}
}
\usage{
ref.grid(object, ...)
    
recover.data(object, ...)
lsm.basis(object, ...)
}

\arguments{
  \item{object}{A model object in a supported class.}
  \item{\dots}{Additional arguments passed to companion functions in the \pkg{emmeans} package.}
} %--- end of \arguments

\value{
\pkg{lsmeans} now passes all its computations to \pkg{emmeans}, and the return values 
are thus what is returned by the corresponding functions \code{\link{ref_grid}}, \code{\link{recover_data}}, and \code{\link{emm_basis}}, respectively.
} 


%%%\references{}

\author{
Russell V. Lenth
}
%\note{}

\examples{
  fiber.lm <- lm(strength ~ machine + diameter, data = fiber)
  rg <- ref.grid(fiber.lm, at = list(diameter = c(20, 24, 28)))
  rg
  
  # Note this is an emmGrid object defined in emmeans. The old "ref.grid"
  # class is now an extension of this:
  r.g. <- new("ref.grid", rg)
  lsmeans(r.g., "machine")
}
% R documentation directory.
\keyword{ regression }
\keyword{ models }% __ONLY ONE__ keyword per line
