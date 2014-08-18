# S3 plot method for lsmobj objects (NOT ref.grid as relies on pri.fac attribute etc.)
# ... are arguments sent to update()


plot.lsmobj = function(x, y, ylab = x@misc$estName, ...) {
    object = update(x, ..., silent = TRUE)
    summ = summary(object, infer = c(TRUE, FALSE))
    
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    
    # Panel functions...
    prepanel.ci = function(x, y, lcl, ucl, subscripts, ...) {
        x = as.numeric(x)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        list(ylim = range(y, ucl, lcl, finite = TRUE)) 
    }
    panel.ci <- function(x, y, lcl, ucl, subscripts, pch = 16, ...) {
        x = as.numeric(x)
        y = as.numeric(y)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        panel.arrows(x, lcl, x, ucl, length = .5, unit = "char",
                     angle = 90, code = 3)
        panel.xyplot(x, y, pch=16, ...)
    }
    
    k = ncol(summ)
    lcl = summ[, k-1]
    ucl = summ[, k]
    
    priv = attr(summ, "pri.vars")
    summ$pri.fac = as.factor(do.call(paste, summ[priv]))
    chform = paste(object@misc$estName, "~", "pri.fac")
    
    byv = attr(summ, "by.vars")
    if (!is.null(byv)) {
        ##summ$by.fac = as.factor(do.call(paste, summ[byv]))
        ##chform = paste(chform, "| by.fac")
        chform = paste(chform, "|", paste(byv, collapse="*"))
    }
    
    form = as.formula(chform)
    
    xyplot(form, prepanel=prepanel.ci, panel=panel.ci, 
           xlab = paste(priv, collapse=":"), ylab = ylab,
           data = summ, lcl=lcl, ucl=ucl, ...)
}