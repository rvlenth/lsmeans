# lsmip code - interaction plots

# object - a model object supported by lsmeans
# formula - a formula of the form  x.factors ~ trace.factors | panel.factors
lsmip = function(object, formula, pch=c(1,2,6,7,9,10,15:20), lty=1, col=NULL, ...) {
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    if (length(formula) < 3)
        stop("'formula' must be two-sided, e.g. trace.factor ~ x.factor")
    ylab = "Least-squares mean"
    
# Glean the parts of ... to use in lsmeans call
    # arguments allowed to be passed
    lsa.allowed = c("at","trend","cov.reduce","fac.reduce")
    xargs = list(...)
    lsmopts = list()
    for (arg in names(xargs)) {
        idx = pmatch(arg, lsa.allowed)
        if (!is.na(idx)) {
            lsmopts[[lsa.allowed[idx]]] = xargs[[arg]]
        }
    }
    
# Get the lsmeans - either from an existing lsm object...    
    if (inherits(object,"lsm"))
        lsms = object[[1]]
    
# ... or from a model        
    else {
        allvars = all.vars(formula)
        lsmopts$object = object
        lsmopts$specs = as.formula(paste("~", paste(allvars, collapse = "+")))
        lsms = do.call("lsmeans", lsmopts) [[1]]
        ylab = paste(ylab, "of", formula(object)[[2]])
    }
    
    
    tvars = all.vars(formula[[2]])
    tv = do.call(paste, lsms[tvars])
    lsms$tvar = factor(tv, levels=unique(tv))
    rhs = strsplit(as.character(formula[3]), "\\|")[[1]]
    xvars = all.vars(as.formula(paste("~", rhs[[1]])))
    xv = do.call(paste, lsms[xvars])
    lsms$xvar = factor(xv, levels = unique(xv))
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
    grobj = xyplot(plotform, groups=~tvar, data=lsms, 
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
