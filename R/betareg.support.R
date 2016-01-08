# Support for 'betareg' class

# NEED TO FIX THIS for precision-related modes
recover.data.betareg = function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)), object$na.action, ...)
}

# PRELIMINARY...
# Currently works correctly only for "resp" and "link" modes
lsm.basis.betareg = function(object, trms, xlev, grid, 
        mode = c("response", "link", "precision", "variance", "quantile"), 
        quantile = .5, ...) {
    mode = match.arg(mode)
    if (mode == "quantile")
        stop('"quantile" mode is not yet supported.')
    
    model = ifelse(mode %in% c("response", "link"), "mean", "precision")
    
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts$mean)
    bhat = coef(object, model = model) 
    V = vcov(object, model = model)
    nbasis = estimability::all.estble
    misc = list(tran = object$link[[model]]$name)
    dffun = function(k, dfargs) NA
    dfargs = list()
    if (mode %in% c("response", "variance")) {
        misc$postGridHook = ".betareg.pg"
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}

# Post-grid hook for simple back-transforming
.betareg.pg = function(object) {
    object@misc$postGridHook = NULL
    regrid(object, transform = TRUE)
}
