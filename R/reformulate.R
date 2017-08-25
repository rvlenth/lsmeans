######################################################################
### Contributed by Jonathon Love, https://github.com/jonathon-love ###
######################################################################

# reformulate for us internally in lsmeans
# same as stats::reformulate, except it surrounds term labels with backsticks

# RVL note: I renamed it .reformulate to avoid certain issues.
#   For example I need reformulate() sometimes to strip off function calls
#   and this .reformulate works quite differently.

.reformulate <- function (termlabels, response = NULL, intercept = TRUE)
{
    if (!is.character(termlabels) || !length(termlabels))
        stop("'termlabels' must be a character vector of length at least one")
    has.resp <- !is.null(response)
    termtext <- paste(if (has.resp)
        "response", "~", paste0("`", termlabels, "`", collapse = "+"), collapse = "")
    if (!intercept)
        termtext <- paste(termtext, "- 1")
    rval <- eval(parse(text = termtext, keep.source = FALSE)[[1L]])
    if (has.resp)
        rval[[2L]] <- if (is.character(response))
            as.symbol(response)
        else response
    environment(rval) <- parent.frame()
    rval
}
