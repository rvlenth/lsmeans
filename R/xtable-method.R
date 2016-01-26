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
        attr(xList, "message") = attr(x, "mesg")
        xList = xtable::xtableList(xList, caption = caption, label = label, 
                            align = align, digits = digits, display = display, 
                            auto = auto, ...)
    }
    else {
        xList = as.data.frame(x)
        xList = xtable::xtable(xList, caption = caption, 
                                   label = label, align = align, digits = digits, display = display, 
                                   auto = auto, ...)
    }
    return(xList)
}