# lsmip code - interaction plots

# object - a model object supported by lsmeans
# formula - a formula of the form  x.factors ~ trace.factors | panel.factors
lsmip = function(object, formula, se.bars = TRUE, ...) {
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    if (length(formula) < 3)
        stop("'formula' must be two-sided, e.g. trace.factor ~ x.factor")
    allvars = all.vars(formula)
    form = as.formula(paste("~", paste(allvars, collapse = "+")))
    lsms = lsmeans(object, form)[[1]]
    tvars = all.vars(formula[[2]])
    lsms$tvar = factor(do.call(paste, lsms[tvars]))
    rhs = strsplit(as.character(formula[3]), "\\|")[[1]]
    xvars = all.vars(as.formula(paste("~", rhs[[1]])))
    lsms$xvar = factor(do.call(paste, lsms[xvars]))
    lsms = lsms[order(lsms$xvar), ]
    plotform = lsmean ~ xvar
    if (length(rhs) > 1) {
        byvars = all.vars(as.formula(paste("~", rhs[[2]])))
        plotform = as.formula(paste("lsmean ~ xvar |", paste(byvars, collapse="*")))
    }
    my.key = function(tvars) 
        list(space="right", 
             title = paste(tvars, collapse=" * "), 
             cex.title=1)
    # The strips the way I want them
    my.strip = function(...)
        strip.default(..., strip.names = c(TRUE,TRUE), sep = " = ")
    # My panel function
    my.panel = function(x, y, subscripts, ucl, lwd, cex, ngrp, group.number, ...) {
        ucl = ucl[subscripts]
        lcl = 2 * y - ucl  # = y - (ucl - y)
        col.line = list(...)$col.line
        offset = .00*(group.number - (ngrp+1)/2)
        panel.xyplot(x, y, subscripts=subscripts, lwd=2, cex=1.5, ...)
        panel.arrows(x+offset, lcl, x+offset, ucl, angle=90, 
                     subscripts=subscripts, code=3, col=col.line, alpha=.10, lwd=10)
    }
    my.main.panel = function(..., groups)
        panel.superpose(..., panel.groups=my.panel, groups=groups,
                        type="o", ngrp = length(unique(groups)))
    my.prepanel = function(x, y, ucl, subscripts, ...) {        
        xlim = range(as.numeric(x))
        del = diff(xlim)/20
        list(ylim = range(c(ucl[subscripts], 2*y-ucl[subscripts]), finite=TRUE))
    }
    grobj = xyplot(plotform, groups= tvar, data=lsms, ucl = lsms$upper.CL,
        xlab = paste("Levels of", paste(xvars, collapse=" * ")),
                   ylab = paste("Least-squares mean of", formula(object)[[2]]),
        strip = my.strip,
        auto.key = my.key(tvars),
        panel = my.main.panel, prepanel = my.prepanel
    )
    print(grobj)
    invisible(lsms)
}