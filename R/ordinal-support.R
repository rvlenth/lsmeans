#--------------------------------------------------------------
### support for the ordinal package
# (initial setup thanks to Maxime Herve, RVAideMemoire package)

recover.data.clm = function(object, ...) {
    if (is.null(object$S.terms) && is.null(object$nom.terms))
        recover.data.lm(object, ...)
    else { # bring-in predictors from loc, scale, and nom models
        trms = delete.response(object$terms)
        #preds = union(all.vars(trms), union(all.vars(object$S.terms), all.vars(object$nom.terms)))
        x.preds = union(all.vars(object$S.terms), all.vars(object$nom.terms))
        #x.trms = update(trms, reformulate(preds))
        x.trms = terms(update(trms, reformulate(c(".", x.preds))))
        recover.data(object$call, x.trms, object$na.action, ...)
    }
}

# For now at least, clmm doesn't cover scale, nominal options
recover.data.clmm = recover.data.lm

# Note: For ALL thresholds, object$Theta has all the threshold values
# for the different cuts (same as object$alpha when threshold=="flexible")
# and object$tJac is s.t. tJac %*% alpha = Theta
# Note also that some functions of cut are constrained to be zero when
# threshold != "flexible". Can get basis using nonest.basis(t(tJac))

lsm.basis.clm = function (object, trms, xlev, grid, ...) {
    contrasts = object$contrasts
    # Remember trms is trumped-up to include scale and nominal predictors.
    # Recover the actual terms for the principal model
    trms = delete.response(object$terms)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    xint = match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) {
        X = X[, -xint, drop = FALSE]
    }
    V = vcov(object)
    k = length(object$alpha)
    # TO DO: Figure out what to do if cuts are aliased...
    # Maybe we just omit certain columns of tJac, alter anm accordingly?
    tJac = object$tJac
    anm = names(object$alpha)
    bnm = names(object$beta)[!is.na(object$beta)]
    cnm = dimnames(tJac)[[1]]
    zeta = object$zeta
    if (!is.null(zeta)) {
        # Overcome a problem whereby some dimensions of V could have the same name!
        names(zeta) = paste(".S",names(zeta), sep=".")
        zkeep = which(!is.na(zeta))
        kz = length(zeta[!is.na(zeta)])
        nm = dimnames(V)[[1]]
        nm[length(nm) - kz + seq_len(kz)] = names(zeta[zkeep])
        dimnames(V) = list(nm, nm)
        bnm = c(bnm, names(zeta[zkeep]))
    }
    allnm = c(cnm, bnm)
    if (is.null(cnm)) {
        cnm = paste(seq_len(nrow(tJac)), "|", 1 + seq_len(nrow(tJac)), sep = "")
    }
    if (object$threshold == "flexible") {
        bhat = c(object$alpha, object$beta)
        V = V[allnm, allnm]
    }
    else {   # We need to transform the alpha part of the coefficients
        bhat = c(tJac %*% object$alpha, object$beta)
        names(bhat)[seq_along(cnm)] = cnm
        V = rbind(tJac %*% V[anm, ], V[bnm, ])
        V = cbind(V[, anm, drop = FALSE] %*% t(tJac), V[, bnm])
        dimnames(V) = list(allnm, allnm)
        k = nrow(tJac)
    }
    j = matrix(1, nrow = k, ncol = 1)
    J = matrix(1, nrow = nrow(X), ncol = 1)
    X = cbind(kronecker(diag(1, k), J), kronecker(-j, X))
    link = as.character(object$info$link)
    misc = list(ylevs = list(cut = cnm), tran = link, 
                inv.lbl = "cumprob")
    if (sum(is.na(object$beta)) > 0) {
        nbasis = nonest.basis(model.matrix(object)$X)
        # Expand the intercept part of this across the cut points
        nbasis = apply(nbasis, 2, function(x) c(-rep(x[1], k), x[-1]))
    }
    else
        nbasis = matrix(NA)
    if (any(object$aliased$alpha))
        message("Note: Estimability is not checked for thresholds")
    
    # deal with scale part of the model
    # TODO: deal with any offsets
    if (!is.null(zeta)) {
        ms = model.frame(object$S.terms, grid, na.action = na.pass, xlev = object$S.xlevels)
        S = model.matrix(object$S.terms, ms, contrasts.arg = object$S.contrasts)
        # keep the right cols (based on unaltered names!)
        S = S[, names(object$zeta), drop = FALSE]
        if (any(is.na(zeta))) {
            nb = nonest.basis(model.matrix(object)$S)
            estble = apply(S, 1, .is.estble, nb)
            S[!estble, ] = NA
            S = S[, zkeep, drop = FALSE]
            zeta = zeta[zkeep]
        }
        offset = 0
        if (!is.null(off.idx <- attr(object$S.terms, "offset"))) {
            offset = rep(0, nrow(grid))
            tvars = attr(object$S.terms, "variables")
            for (i in off.idx)
                offset = offset + eval(tvars[[i+1]], grid)
        }
        sigma = exp(offset + S %*% zeta)
        D = diag(rep(1/sigma, k))
        X = D %*% X
        # OK, now to get the SEs right, we need to transform to the X%*%beta space
        # 1. add -S * eta to X
        bkeep = which(!is.na(bhat))
        eta = as.numeric(X[, bkeep] %*% bhat[bkeep])
        offset = rep(0, nrow(grid))
        if (!is.null(off.idx <- attr(trms, "offset"))) {
            tvars = attr(trms, "variables")
            for (i in off.idx)
                offset = offset + eval(tvars[[i+1]], grid)
        }
        doff = as.numeric(D %*% rep(offset, k))
        L = cbind(X, diag(eta + doff) %*% kronecker(matrix(-1,nrow=k), S))
        V = L %*% V %*% t(L)
        X = diag(1, nrow(L))
        bhat = eta - doff # compensate for fact that offset value is still gonna be in grid
        # TODO: Can we figure out estimability?        
    }
    
    dffun = function(...) NA
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}

lsm.basis.clmm = function (object, trms, xlev, grid, ...) {
    if(is.null(object$Hessian)) {
        message("Updating the model to obtain the Hessian...")
        object = update(object, Hess = TRUE)
    }
    # borrowed from Maxime's code -- need to understand this better, e.g. when it happens
    H = object$Hessian
    if (any(apply(object$Hessian, 1, function(x) all(x == 0)))) {
        H = H[names(coef(object)), names(coef(object))]
        object$Hessian = H
    }
    lsm.basis.clm(object, trms, xlev, grid, ...)
}
