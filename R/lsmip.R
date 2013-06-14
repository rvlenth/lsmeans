# lsmip code - interaction plots

# object - a model object supported by lsmeans
# formula - a formula of the form  x.factors ~ trace.factors | panel.factors
lsmip = function(object, formula, pch=c(1,2,6,7,9,10,15:20), lty=1, col=NULL, ...) {
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    if (length(formula) < 3)
        stop("'formula' must be two-sided, e.g. trace.factor ~ x.factor")
    ylab = "Least-squares mean"
    if (inherits(object,"lsm"))
        lsms = object[[1]]
    else {
        allvars = all.vars(formula)
        form = as.formula(paste("~", paste(allvars, collapse = "+")))
        lsms = lsmeans(object, form)[[1]]
        ylab = paste(ylab, "of", formula(object)[[2]])
    }
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
             points = TRUE, 
             lines=length(lty) > 1,
             cex.title=1)
    # The strips the way I want them
    my.strip = function(...)
        strip.default(..., strip.names = c(TRUE,TRUE), sep = " = ")
    TP = TP.orig = trellis.par.get()
    TP$superpose.symbol$pch = pch
    TP$superpose.line$lty = lty
    if (!is.null(col)) TP$superpose.symbol$col = TP$superpose.line$col = col
    trellis.par.set(TP)
    grobj = xyplot(plotform, groups= tvar, data=lsms, 
                   xlab = paste("Levels of", paste(xvars, collapse=" * ")),
                   ylab = ylab,
                   strip = my.strip,
                   auto.key = my.key(tvars), 
                   type=c("p","l"), 
                   ... )
    print(grobj)
    trellis.par.set(TP.orig)
    invisible(lsms)
}