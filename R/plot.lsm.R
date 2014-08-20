# S3 plot method for lsmobj objects (NOT ref.grid as relies on pri.fac attribute etc.)
# ... are arguments sent to update()


plot.lsmobj = function(x, y, type, ...) {
    if(!missing(type))
        object = update(x, predict.type = type, ..., silent = TRUE)
    else
        object = update(x, ..., silent = TRUE)
    summ = summary(object, infer = c(TRUE, FALSE))
    estName = attr(summ, "estName")
    plot(summ, ...)
}

# May use in place of plot.lsmobj but no control over level etc.

plot.summary.ref.grid = function(x, y, horizontal = TRUE, xlab, ylab, ...) {
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    
    summ = x # so I don't get confused
    estName = attr(summ, "estName")
    clNames = attr(summ, "clNames")
    if (is.null(clNames)) {
        warning("No information available to display confidence limits")
        lcl = ucl = summ[[estName]]
    }
    else {
        lcl = summ[[clNames[1]]]
        ucl = summ[[clNames[2]]]
    }
    
    # Panel functions...
    prepanel.ci = function(x, y, horizontal=TRUE, lcl, ucl, subscripts, ...) {
        x = as.numeric(x)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        if (horizontal)
            list(xlim = range(x, ucl, lcl, finite = TRUE)) 
        else
            list(ylim = range(y, ucl, lcl, finite = TRUE)) 
    }
    panel.ci <- function(x, y, horizontal=TRUE, lcl, ucl, subscripts, pch = 16, 
                         lty = dot.line$lty, lwd = dot.line$lwd, 
                         col = dot.symbol$col, col.line = dot.line$col, ...) {
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        x = as.numeric(x)
        y = as.numeric(y)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        if(horizontal) {
            panel.abline(h = unique(y), col = col.line, lty = lty, lwd = lwd)
            panel.arrows(lcl, y, ucl, y, col=col, length = .5, unit = "char", angle = 90, code = 3)
        }
        else {
            panel.abline(v = unique(x), col = col.line, lty = lty, lwd = lwd)
            panel.arrows(x, lcl, x, ucl, col=col, length = .5, unit = "char", angle = 90, code = 3)
        }
        panel.xyplot(x, y, pch=16, ...)
    }
    my.strip = strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
    
    priv = attr(summ, "pri.vars")
    pf = do.call(paste, summ[priv])
    summ$pri.fac = factor(pf, levels=unique(pf))
    chform = ifelse(horizontal,
                    paste("pri.fac ~", estName),
                    paste(estName, "~ pri.fac"))
    
    byv = attr(summ, "by.vars")
    if (!is.null(byv)) {
        chform = paste(chform, "|", paste(byv, collapse="*"))
    }
    
    facName = paste(priv, collapse=":")
    form = as.formula(chform)
    if (horizontal) {
        if (missing(xlab)) xlab = estName
        if (missing(ylab)) ylab = facName
        dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = TRUE,
                ylab = ylab, xlab = xlab,
                data = summ, lcl=lcl, ucl=ucl, ...)
    }
    else {
        if (missing(xlab)) xlab = facName
        if (missing(ylab)) ylab = estName
        dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = FALSE,
                xlab = paste(priv, collapse=":"), ylab = ylab,
                data = summ, lcl=lcl, ucl=ucl, ...)
    }
    
    
}
