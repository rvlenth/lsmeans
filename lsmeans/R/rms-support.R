# Support for objects in the *rms* package

recover.data.rms = function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)), object$na.action$omit, ...)
}

# TODO: 
# 1. If multivariate - like mlm method?
# 2. orm cases?

lsm.basis.rms = function(object, trms, xlev, grid, mode = c("middle", "latent", "linear.predictor"), ...) {
    mode = match.arg(mode)
    bhat = coef(object) 
    V = vcov(object, intercepts = "all")
    misc = list()
    
    X = predict(object, newdata = grid, type = "x")
    #xnames = dimnames(X)[[2]]
    #intcpts = setdiff(names(bhat), xnames)
    nint = length(bhat) - ncol(X)
    intcpts = names(bhat)[seq_len(nint)]
    xnames = setdiff(names(bhat), intcpts)
        
    if (length(intcpts) == 1) 
        mode = "single" # stealth mode for ordinary single-intercept case
    if (mode %in% c("single", "middle", "latent")) {
        X = cbind(1, X)
        mididx = ifelse(mode != "middle", 1, as.integer((1 + length(intcpts)) / 2))
        dimnames(X)[[2]][1] = switch(mode,
            single = intcpts,
            middle = intcpts[mididx],
            latent = "avg.intercept")
        if (mode == "middle") {
            nms = c(intcpts[mididx], xnames)
            bhat = bhat[nms]
            V = V[nms, nms, drop = FALSE]
        }
        else if (mode == "latent") {
            bhat = c(mean(bhat[intcpts]), bhat[xnames])
            nx = length(xnames)
            J1 = rbind(rep(1/nint, nint), 
                       matrix(0, nrow = nx, ncol = nint))
            J2 = rbind(0, diag(1, nx))
            J = cbind(J1, J2)
            V = J %*% V %*% t(J)
        }
        ### else mode == "single" and all is OK as it is
    }
    else { # mode == "linear.predictor"
        misc$ylevs = list(cut = intcpts)
        I = diag(1, nint)
        J = matrix(1, nrow = nrow(X))
        JJ = matrix(1, nrow=nint)
        X = cbind(kronecker(I, J), kronecker(JJ, X))
        # Note V is correct as-is
        dimnames(X)[[2]] = c(intcpts, xnames)
    }
    
    # I think rms does not allow rank deficiency...
    nbasis = estimability::all.estble
    if (inherits(object, "glm")) {
        misc = .std.link.labels(object$family, misc)
        dffun = function(k, dfargs) NA
        dfargs = list()
    }
    else {
        dfargs = list(df = object$df.residual)
        if (is.null(dfargs$df)) 
            dfargs$df = NA
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, 
         dffun=dffun, dfargs=dfargs, misc=misc)
}
