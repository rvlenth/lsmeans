\name{glht.lsmobj}
\alias{glht.lsmobj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
\pkg{lsmeans} support for \code{glht}
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glht.lsmobj(model, linfct, by, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{linfct}{
%%     ~~Describe \code{linfct} here~~
}
  \item{by}{
%%     ~~Describe \code{by} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, linfct, by, ...) 
{
    object = linfct
    args = list(model = model, ...)
    if (is.null(args$df)) {
        df = summary(linfct)$df
        if (any(!is.na(df))) {
            args$df = max(1, as.integer(mean(df, na.rm = TRUE) + 
                0.25))
            message("Note: df set to ", args$df)
        }
    }
    if (missing(by)) 
        by = object@misc$by.vars
    if (is.null(by)) {
        args$linfct = object@linfct
        return(do.call("glht", args))
    }
    by.rows = .find.by.rows(object@grid, by)
    result = lapply(by.rows, function(r) {
        args$linfct = object@linfct[r, , drop = FALSE]
        do.call("glht", args)
    })
    bylevs = lapply(by, function(byv) unique(object@grid[[byv]]))
    names(bylevs) = by
    bygrid = do.call("expand.grid", bylevs)
    levlbls = lapply(by, function(byv) paste(byv, "=", bygrid[[byv]]))
    levlbls$sep = ", "
    names(result) = do.call("paste", levlbls)
    class(result) = c("glht.list", "list")
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
