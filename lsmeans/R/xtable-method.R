### xtable method
# Modified from xtableLSMeans function provided by David Scott

xtable.ref.grid = function(x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, auto = FALSE, ...) 
{
    xtable.summary.ref.grid(summary(x, ...), caption = caption, label = label, align = align, digits = digits, 
           display = display, auto = auto)
}

xtable.summary.ref.grid = function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
          display = NULL, auto = FALSE, ...) 
{
    if (!is.null(byv <- attr(x, "by.vars"))) {
        byc = which(names(x) %in% byv)
        xList = split(as.data.frame(x), f = x[, byc])
        labs = rep("", length(xList))
        for (i in 1:length(xList)) {
            levs = sapply(xList[[i]][1, byc], as.character)
            labs[i] = paste(paste(byv, levs, sep = " = "), collapse = ", ")
            xList[[i]] = as.data.frame(xList[[i]][, -byc, drop = FALSE])
        }
        attr(xList, "subheadings") = labs
    }
    else {
        xList = list(as.data.frame(x))
    }
    attr(xList, "message") = attr(x, "mesg")
    result = xtable::xtableList(xList, caption = caption, label = label, 
       align = align, digits = digits, display = display, 
       auto = auto, ...)
    class(result) = c("xtable.lsm", "xtableList")
    result
}

# My own print method
print.xtable.lsm = function(x, include.rownames = FALSE, 
                            sanitize.message.function = footnotesize,
                            ...)
{
    footnotesize = function(x) paste0("{\\footnotesize ", x, "}")
    invisible(print.xtableList(x, include.rownames = include.rownames, sanitize.message.function = sanitize.message.function, ...))
}    