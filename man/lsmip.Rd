\name{lsmip}
\alias{lsmip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-squares means interaction plot
}
\description{
This function creates an interaction plot of the least-squares means based on a fitted model and a simple formula specification.
}
\usage{
lsmip(object, formula, pch = c(1,2,6,7,9,10,15:20), lty = 1, col = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{lsm}, or a fitted model of a class supported by \code{\link{lsmeans}}.
}
  \item{formula}{
Formula of the form \code{trace.factors ~ x.factors | by.factors}. The least-squares means are plotted against \code{x.factor} for each level of \code{trace.factors}. \code{by.factors} is optional, but if present, it determines separate panels. Each element of this formula may be a single factor in the model, or a combination of factors using the \code{*} operator.
}
  \item{pch}{
The plotting characters to use for each group (i.e., levels of \code{trace.factors}). They are recycled as needed.
}
  \item{lty}{
The line types to use for each group. Recycled as needed.
}
  \item{col}{
The colors to use for each group, recycled as needed. If not specified,
the default trellis colors are used.
}
  \item{\dots}{
Additional arguments passed to \code{\link{lsmeans}} or to \code{\link[lattice]{xyplot}}. 
}
}
\details{
If \code{object} is a fitted model, \code{\link{lsmeans}} is called with an appropriate specification to obtain least-squares means for each combination of the factors present in \code{formula} (in addition, any arguments in \code{\dots} that match \code{at}, \code{trend}, \code{cov.reduce}, or \code{fac.reduce} are passed to \code{lsmeans}).
Otherwise, if \code{object} is an \code{lsm} object, its first element is used, and it must contain one \code{lsmean} value for each combination of the factors present in \code{formula}.
}
\value{
(Invisibly), the table of least-squares means that were plotted.
}
\author{
Russell V. Lenth, The University of Iowa
}
\note{
This function uses the \code{\link[lattice]{xyplot}} function in the \code{lattice} package (an error is returned if \code{lattice} is not installed). Conceptually, it is equivalent to \code{\link{interaction.plot}} where the summarization function is the least-squares means.
}

\seealso{
\code{\link{interaction.plot}}
}
\examples{
#--- Two-factor example
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)

# Following plot is the same as the usual interaction plot of the data
lsmip(warp.lm, wool ~ tension)

#--- Four-factor example from Box, Hunter, Hunter (2nd ed) p. 199
process <- expand.grid(Catalyst=factor(c(10,15)), Temp=factor(c(220,240)), 
                       Pressure=factor(c(50,80)), Conc=factor(c(10,12)))
process$Conversion <- c(70,60,89,81,69,62,88,81,60,49,88,82,60,52,86,79)
process.lm <- lm(Conversion ~ Catalyst*Temp*Pressure*Conc, data = process)

op <- options(warn=-1) # zero df for error, suppress warnings

# Multi-panel interaction plot
lsmip(process.lm, Catalyst ~ Temp | Pressure*Conc)

# One panel using different factor combinations and custom symbols and lines 
lsmip(process.lm, Catalyst*Temp ~ Pressure*Conc, 
      lty=1:2, col=1:2, pch=c(1,1,2,2))

# Simplified model with just one two-way interaction
process.lm2 <- update(process.lm, . ~ Catalyst + Pressure + Temp*Conc)
lsmip(process.lm2, Catalyst*Temp ~ Pressure*Conc, 
      lty=1:2, col=1:2, pch=c(1,1,2,2))

options(op)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ regression }
