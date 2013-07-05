\name{cld.lsm}
\alias{cld}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cld
}
\description{
  cld provides the results of lsmeans with a compact display letters for the
  differences between the means or trends.
}
\usage{
cld(obj, reversed = FALSE, omitNotSig  = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
    a object of class lsm
}
  \item{revered}{
    reverse the letters order, i.e (cba) instead of (abc)
}
  \item{omitNotSig}{
    omit not signifactive letters, i.e. if there are no differences inside 
    a group then the letters for that group become "" instead of all "a".
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  It returns a data.frame equal to results of lsmeans plus the cld.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 Luciano Selzer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
  cld(lsmeans(warp.lm,  pairwise ~ tension|wool))
}
