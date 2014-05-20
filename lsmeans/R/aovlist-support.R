# lsmeans support for aovlist objects

recover.data.aovlist = function(object, ...) {
    fcall = attr(object, "call")
    trms = terms(object)
    # Find the Error terms
    lbls = attr(trms, "term.labels")
    err.idx = grep("^Error\\(", lbls)
    newf = as.formula(paste(c(".~.", lbls[err.idx]), collapse = "-"))
    trms = terms(update(trms, newf))
    recover.data(fcall, delete.response(trms), na.action = attr(object, "na.action"), ...)
}

# This works great for balanced experiments, and goes horribly wrong
# even for slightly unbalanced ones. So I abort on these kinds of cases
lsm.basis.aovlist = function (object, trms, xlev, grid) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr = attr(object, "contrasts")
    X = model.matrix(trms, m, contrasts.arg = contr)
    xnms = dimnames(X)[[2]]
    
    # Check for situations we can't handle...
    nms = unlist(lapply(object, function(x) names(coef(x))))
    if(max(table(match(nms, xnms))) > 1)
        stop("Can't handle situations where the same effects appear in more than one stratum")
    message("NOTE: 'aovlist' support is experimental.")
    if(!all(unlist(contr) %in% c("contr.helmert","contr.poly")))
        warning("Model uses non-orthogonal contrasts - lsmeans results probably incorrect\n",
            "Contrast results probably OK if on one stratum (integer df)")
    
    # initialize arrays
    k = length(xnms)
    bhat = rep(0, k)
    V = matrix(0, nrow=k, ncol=k)
    names(bhat) = xnms
    dimnames(V) = list(xnms, xnms)
    Vmats = Vidx = Vdf = list()
    for (nm in setdiff(names(object), "(Intercept)")) {
        x = object[[nm]]
        bi = coef(x)
        bi = bi[!is.na(bi)]
        ii = Vidx[[nm]] = match(names(bi), xnms)
        bhat[ii] = bi
        Vmats[[nm]] = vcov(x)
        Vdf[[nm]] = x$df
        V[ii,ii] = Vmats[[nm]]
    }
    x <- object[["(Intercept)"]]
    if (!is.null(x)) {
        # The intercept belongs in the 1st error stratum
        # So have to add a row and column to its covariance matrix
        bhat[1] = x$coefficients[1]
        Vidx[[1]] = ii = c(1, Vidx[[1]])
        k = length(ii)
        vv = matrix(0, nrow=k, ncol=k)
        if (k > 1) vv[2:k,2:k] = Vmats[[1]]
        # Variance of intercept is EMS of this stratum divided by N
        # Here I'm assuming there are no weights
        N = sum(sapply(object, function(x) length(x$residuals)))
        V[1,1] = vv[1,1] = sum(object[[2]]$residuals^2) / object[[2]]$df / N
        #dimnames(vv) = list(c(xnms[ii], xnms[ii]))
        Vmats[[1]] = vv
    }
    dfargs = list(Vmats=Vmats, Vidx=Vidx, Vdf=unlist(Vdf))
    dffun = function(k, dfargs) {
        v = sapply(1:length(dfargs$Vdf), function(j) {
            ii = dfargs$Vidx[[j]]
            sum(k[ii] * dfargs$Vmats[[j]] %*% k[ii])
        })
        sum(v)^2 / sum(v^2 / dfargs$Vdf) # Good ole Satterthwaite
    }
    nbasis = matrix(NA)  # Consider this further?
    misc = list()
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs, misc = misc)
}