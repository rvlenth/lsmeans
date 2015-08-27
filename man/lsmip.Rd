\name{lsmip}
\alias{lsmip}
\alias{lsmip.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Least-squares means interaction plot
}
\description{
This function creates an interaction plot of the least-squares means based on a fitted model and a simple formula specification.
}
\usage{
\method{lsmip}{default}(object, formula, type,
    pch = c(1,2,6,7,9,10,15:20), 
    lty = 1, col = NULL, plotit = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{lsmobj}, or a fitted model of a class supported by \code{\link{lsmeans}}.
}
  \item{formula}{
Formula of the form \code{trace.factors ~ x.factors | by.factors}. The least-squares means are plotted against \code{x.factor} for each level of \code{trace.factors}. \code{by.factors} is optional, but if present, it determines separate panels. Each element of this formula may be a single factor in the model, or a combination of factors using the \code{*} operator.
}
\item{type}{
As in \code{\link[=predict.ref.grid]{predict}}, this determines whether we want to inverse-transform the predictions (\samp{type="response"}) or not (any other choice). The default is \code{"link"}, unless the \code{"predict.type"} option is in force; see \code{\link{lsm.options}}.
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
  \item{plotit}{
If \code{TRUE}, the plot is displayed. Otherwise, one may use the \code{"lattice"} attribute of the returned object and print it, perhaps after additional manipulation.  
}
  \item{\dots}{
Additional arguments passed to \code{\link{lsmeans}} or to \code{\link[lattice]{xyplot}}. 
}
}
\details{
If \code{object} is a fitted model, \code{\link{lsmeans}} is called with an appropriate specification to obtain least-squares means for each combination of the factors present in \code{formula} (in addition, any arguments in \code{\dots} that match \code{at}, \code{trend}, \code{cov.reduce}, or \code{fac.reduce} are passed to \code{lsmeans}).
Otherwise, if \code{object} is an \code{lsmobj} object, its first element is used, and it must contain one \code{lsmean} value for each combination of the factors present in \code{formula}.
}
\value{
(Invisibly), a \code{\link{data.frame}} with the table of least-squares means that were plotted, with an additional \code{"lattice"} attribute containing the \code{trellis} object for the plot.
}
\author{
Russell V. Lenth
}
\note{
This function uses the \code{\link[lattice]{xyplot}} function in the \code{lattice} package (an error is returned if \code{lattice} is not installed). Conceptually, it is equivalent to \code{\link{interaction.plot}} where the summarization function is the least-squares means.
}

\seealso{
\code{\link{interaction.plot}}
}
\examples{
require(lsmeans)
require(lattice)

#--- Two-factor example
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)

# Following plot is the same as the usual interaction plot of the data
lsmip(warp.lm, wool ~ tension)

#--- Three-factor example
noise.lm = lm(noise ~ size * type * side, data = auto.noise)

# Separate interaction plots of size by type, for each side
lsmip(noise.lm, type ~ size | side)

# One interaction plot, using combinations of size and side as the x factor
lsmip(noise.lm, type ~ side * size)

# One interaction plot using combinations of type and side as the trace factor
# customize the colors, line types, and symbols to suggest these combinations
lsmip(noise.lm, type * side ~ size, lty=1:2, col=1:2, pch=c(1,1,2,2))

# 3-way interaction is significant, but doesn't make a lot of visual difference...
noise.lm2 = update(noise.lm, . ~ . - size:type:side)
lsmip(noise.lm2, type * side ~ size, lty=1:2, col=1:2, pch=c(1,1,2,2))
}

\keyword{ models }
\keyword{ regression }
