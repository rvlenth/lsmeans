%    Copyright (c) 2012-2016 Russell V. Lenth                                %

\name{transition}
\alias{transition}
\alias{emmeans-transition}

\title{Transition to \pkg{emmeans}}

\description{
The \pkg{lsmeans} package is being deprecated and further development will take
place in its successor, \pkg{emmeans}. Users may use \pkg{emmeans} in almost 
exactly the same way as \pkg{lsmeans}, but a few function names and internal 
details are changed.}


\details{
In transitioning to \pkg{emmeans}, users will find that the vignettes are
constructed quite differently and that, in those and in the documemntation,
emphasis is placed on \dQuote{estimated marginal means} rather than 
\dQuote{least-squares means}. The latter term is broader and more appropriate 
for use with some models, e.g. ordinal regression, that don't really involve
least-squares methods. That is the reason for the change.

Accordingly, \pkg{emmeans} users are encouraged to use the functions \code{emmeans()}, 
\code{emtrends()}, \code{emmip()}, etc. in lieu of \code{lsmeans()}, etc.
The latter functions \emph{are still available} in \pkg{emmeans}; they run
the corresponding \code{emmxxxx} function and relabel the results.

The \pkg{emmeans} package provides some functions that help convert scripts and R Markdown
files containing \pkg{lsmeans} code so they will work in \pkg{emmeans}.
There is also a function to convert \code{ref.grid} and \code{lsmobj} objects
to the \code{emmGrid} objects used in \pkg{emmeans}.
More extensive information is given in 
\code{vignette("transition-from-lsmeans", package = "emmeans")}.}


\author{Russell V. Lenth}
