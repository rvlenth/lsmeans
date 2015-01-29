# experimental support for nls objects

recover.data.nls = function(object, ...) {
    fcall = object$call
    trms = terms(reformulate(names(object$dataClasses)))
    recover.data(fcall, trms, object$na.action, ...)
}

lsm.basis.nls = function(object, trms, xlev, grid, ...) {
    Vbeta = vcov(object)
    env = object$m$getEnv()
    for (nm in names(grid)) env[[nm]] = grid[[nm]]
    pars = object$m$getAllPars()
    DD = deriv(object$m$formula(), names(pars))
    ests = eval(DD, env)
    bhat = as.numeric(ests)
    grad = attr(ests, "gradient")
    V = grad %*% Vbeta %*% t(grad)
    X = diag(1, nrow(grid))
    list(X=X, bhat=bhat, nbasis=all.estble, V=V, 
         dffun=function(k, dfargs) NA, dfargs=list(), 
         misc=list())
}
    
