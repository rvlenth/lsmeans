# S3 plot method for lsmobj objects (NOT ref.grid as relies on pri.fac attribute etc.)
# ... are arguments sent to update()


plot.lsmobj = function(x, y, ylab = estName, type, ...) {
    if(!missing(type))
        object = update(x, predict.type = type, ..., silent = TRUE)
    else
        object = update(x, ..., silent = TRUE)
    summ = summary(object, infer = c(TRUE, FALSE))
    estName = attr(summ, "estName")
    plot(summ, ylab=ylab, ...)
}

# May use in place of plot.lsmobj but no control over level etc.
plot.summary.ref.grid = function(x, y, ylab = estName, ...) {
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    
    summ = x # so I don't get confused
    estName = attr(summ, "estName")
    clNames = attr(summ, "clNames")
    if (is.null(clNames)) {
        message("No information available to display confidence limits")
        lcl = ucl = summ[[estName]]
    }
    else {
        lcl = summ[[clNames[1]]]
        ucl = summ[[clNames[2]]]
    }
    
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
    
    priv = attr(summ, "pri.vars")
    pf = do.call(paste, summ[priv])
    summ$pri.fac = factor(pf, levels=unique(pf))
    chform = paste(estName, "~", "pri.fac")
    
    byv = attr(summ, "by.vars")
    if (!is.null(byv)) {
        chform = paste(chform, "|", paste(byv, collapse="*"))
    }
    
    form = as.formula(chform)
    
    xyplot(form, prepanel=prepanel.ci, panel=panel.ci, 
           xlab = paste(priv, collapse=":"), ylab = ylab,
           data = summ, lcl=lcl, ucl=ucl, ...)
}

