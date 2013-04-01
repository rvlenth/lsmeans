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
    grobj = xyplot(plotform, groups=~ tvar, type="o", data=lsms,
        xlab = paste("Levels of", paste(xvars, collapse=" * ")),
                   ylab = paste("Least-squares mean of", formula(object)[[2]]),
        auto.key = list(space="right", 
            title = paste(tvars,collapse=" * "), cex.title=1)
    )
    print(grobj)
    invisible(lsms)
}