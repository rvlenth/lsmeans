# Support for objects in the *rms* package

recover.data.rms = function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)), object$na.action$omit, ...)
}

# TODO: 
# 1. If multivariate - like mlm method?
# 2. orm cases?

lsm.basis.rms = function(object, trms, xlev, grid, ...) {
    X = predict(object, newdata = grid, type = "x")
    X = cbind(Intercept = 1, X)
    bhat = coef(object) 
    V = vcov(object)
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    misc = list()
    if (inherits(object, "glm")) {
        misc = .std.link.labels(object$family, misc)
        dffun = function(k, dfargs) NA
        dfargs = list()
    }
    else {
        dfargs = list(df = object$df.residual)
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}
