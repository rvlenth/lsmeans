# Support for 'betareg' class

# mode is same as 'type' in predict.betareg, PLUS 
# mode = "phi.link" refers to link function before back-transforming to "precision"

recover.data.betareg = function(object, mode = c("response", "link", "precision", "phi.link", "variance", "quantile"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode  %in% c("response", "link"))
        mode = "mean"
    if (mode == "phi.link") 
        mode = "precision"
    if(mode %in% c("response", "link", "precision"))
        trms = delete.response(terms(object, model = mode))
    else
        trms = delete.response(object$terms$full)
    # Make sure there's an offset function available
    env = new.env(parent = attr(trms, ".Environment"))
    env$offset = function(x) x
    attr(trms, ".Environment") = env
    recover.data(fcall, trms, object$na.action, ...)
}

# PRELIMINARY...
# Currently works correctly only for "resp", "link", "precision", "phi" modes
lsm.basis.betareg = function(object, trms, xlev, grid, 
        mode = c("response", "link", "precision", "phi.link", "variance", "quantile"), 
        quantile = .5, ...) {
    mode = match.arg(mode)
    if (mode %in% c("variance", "quantile"))
        stop(paste0('"', mode, '" mode is not yet supported.'))
    
    # figure out which parameters we need
    model = if (mode %in% c("response", "link")) "mean"
        else if (mode %in% c("precision", "phi.link")) "precision"
        else "full"
    V = vcov(object, model = model)
    bhat = coef(object, model = model)
    
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts[[model]])
    
    nbasis = estimability::all.estble
    misc = list(tran = object$link[[model]]$name)
    dffun = function(k, dfargs) NA
    dfargs = list()
    if (mode %in% c("response", "precision")) {
        misc$postGridHook = ".betareg.pg"
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}

# Post-grid hook for simple back-transforming
.betareg.pg = function(object) {
    object@misc$postGridHook = NULL
    regrid(object, transform = TRUE)
}


### predict methods
# link: X%*%beta + off_m
# response: mu = h_m(link)
# 
# phi.link: Z%*%gamma + off_p
# precision: phi = h_p(phi.link)
# 
# variance: mu*(1 - mu) / (1 + phi)
# quantile: qbeta(p, mu*phi, (1 - mu)*phi)
#
# Defns:
#   phi = a + b
#   mu = a / (a + b)
# so that phi*mu = a and phi*(1 - mu) = b,
#   Variance = ab / [(a + b)^2 * (a + b + 1)]
