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
#
# opt arg 'latent' - if true, we work with +x'beta, else we use theta_j - x'beta

# NEW VERSION - assumes everything is est'ble
# TODO: figure out non-est basis

lsm.basis.clm = function (object, trms, xlev, grid, 
                          latent = TRUE, rescale = c(0,1), ...) {
    # general stuff
    bhat = coef(object)
    V = vcov(object)
    tJac = object$tJac
    dffun = function(...) NA
    link = as.character(object$info$link)
    cnm = dimnames(object$tJac)[[1]]
    if (is.null(cnm))
        cnm = paste(seq_len(nrow(tJac)), "|", 1 + seq_len(nrow(tJac)), sep = "")
    misc = list()
    
    
    # My strategy is to piece together the needed matrices for each threshold parameter
    # Then assemble the results
    
    ### ----- Location part ----- ###
    contrasts = object$contrasts
    # Remember trms was trumped-up to include scale and nominal predictors.
    # Recover the actual terms for the principal model
    trms = delete.response(object$terms)
    m = model.frame(trms, grid, na.action = na.pass, xlev = object$xlevelssummarundebui)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    xint = match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) {
        X = X[, -xint, drop = FALSE]
    }
    
    ### ----- Nominal part ----- ###
    if (is.null(object$nom.terms))
        NOM = matrix(1, nrow = nrow(X))
    else {
        mn = model.frame(object$nom.terms, grid, na.action = na.pass, xlev = object$nom.xlevels)
        NOM = model.matrix(object$nom.terms, mn, contrasts.arg = object$nom.contrasts)
    }
    bigNom = kronecker(tJac, NOM)
    # cols are in wrong order... I'll get the indexes by transposing a matrix of subscripts
    if (ncol(NOM) > 1)
        bigNom = bigNom[, as.numeric(t(matrix(seq_len(ncol(bigNom)), nrow=ncol(NOM))))]
    
    ### ----- Scale part ----- ###
    if (!is.null(object$S.terms)) {
        ms = model.frame(object$S.terms, grid, na.action = na.pass, xlev = object$S.xlevels)
        S = model.matrix(object$S.terms, ms, contrasts.arg = object$S.contrasts)
        S = S[, names(object$zeta), drop = FALSE]
        if (!is.null(attr(object$S.terms, "offset"))) {
            soff = .get.offset(object$S.terms, grid)
            # we'll add a column to S and adjust bhat and V accordingly
            S = cbind(S, offset = soff)
            bhat = c(bhat, offset = 1)
            V = rbind(cbind(V, offset = 0), offset = 0)
        }
        si = misc$scale.idx = length(object$alpha) + length(object$beta) + seq_len(ncol(S))
        # Make sure there are no name clashes
        names(bhat)[si] = paste(".S", names(object$zeta), sep=".")
        misc$estHook = as.name(".clm.estHook")
        misc$vcovHook = as.name(".clm.vcovHook")
    }
    else
        S = NULL
    
    ### ----- Get non-estimability basis ----- ###
    nbasis = snbasis = matrix(NA)
    if (any(is.na(bhat))) {
        mm = model.matrix(object)
        if (any(is.na(c(object$alpha, object$beta)))) {
            NOMX = mm$X
            if (is.null(mm$NOM))
                NOMX = mm$X
            else
                NOMX = cbind(mm$NOM, mm$X[, -1, drop=false])
            nbasis = nonest.basis(NOMX)
            # replicate the NOM parts
            nomcols = seq_len(ncol(NOM))
            nbasis = apply(nbasis, 2, function(x)
                c(rep(x[nomcols], each = length(object$alpha)), x[-nomcols]))
        }
        if (!is.null(mm$S)) {
            if (any(is.na(object$zeta))) {
                snbasis = nonest.basis(mm$S)
                # put intercept part at end
                snbasis = rbind(snbasis[-1, , drop=FALSE], snbasis[1, ])
                if (!is.null(attr(object$S.terms, "offset")))
                    snbasis = rbind(snbasis, 0)
                snbasis = rbind(matrix(0, ncol=ncol(snbasis), nrow=min(si)-1), snbasis)
                 # Note scale intercept is included, so tack it on to the end of everything
                S = cbind(S, .S.intcpt = 1)
                bhat = c(bhat, .S.intcpt = 0)
                V = rbind(cbind(V, .S.intcpt = 0), .S.intcpt = 0)
                si = misc$scale.idx = c(si, 1 + max(si))
            }
        }
        if (is.na(nbasis[1])) # then only nonest part is scale
            nbasis = snbasis
        else if (!is.na(snbasis[1])) # then have both
            nbasis = cbind(rbind(nbasis, matrix(0, nrow=length(si), ncol=ncol(nbasis))), snbasis)
        # else nbasis stands as-is
    }
    
    if (latent) {
        # Create constant columns for means of scale and nominal parts
        J = matrix(1, nrow = nrow(X))
        nomm = rescale[2] * apply(bigNom, 2, mean)
        X = rescale[2] * X
        if (!is.null(S)) {
            sm = apply(S, 2, mean)
            X = cbind(X, kronecker(-J, matrix(sm, nrow = 1)))
        }
        bigX = cbind(kronecker(-J, matrix(nomm, nrow = 1)), X)
        misc$offset.mult = misc$offset.mult * rescale[2]
        intcpt = seq_len(ncol(tJac))
        bhat[intcpt] = bhat[intcpt] - rescale[1] / rescale[2]
    }
    else { ### ----- Piece together big matrix for each threshold ----- ###
        misc$ylevs = list(cut = cnm)
        misc$tran = link
        misc$inv.lbl = "cumprob"
        misc$offset.mult = -1
        if (!is.null(S))
            X = cbind(X, S)
        J = matrix(1, nrow=nrow(tJac))
        bigX = cbind(bigNom, kronecker(-J, X))
    }
    
    dimnames(bigX)[[2]] = names(bhat)
    
    list(X = bigX, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}

# replacement estimation routine for cases with a scale param
.clm.estHook = function(object, do.se = TRUE, tol = 1e-8) {
    scols = object@misc$scale.idx
    linfct = object@linfct
    bhat = object@bhat
    active = !is.na(bhat)
    bhat[!active] = 0
    if (is.null(tol)) 
        tol = 1e-8
    
    result = sapply(seq_len(nrow(linfct)), function(i) {
        x = linfct[i, ]
        if (.is.estble(x, object@nbasis, tol)) {
            rsigma = exp(sum(linfct[i, scols] * bhat[scols]))
            est = rsigma * sum(bhat[-scols] * x[-scols])
            if (!is.null(object@grid$.offset.))
                est = est + rsigma * object@grid$.offset.[i] # offset is already negated via offset.mult
            if(do.se) {
                x[scols] = est * x[scols] / rsigma
                se = sqrt(.qf.non0(object@V, rsigma * x[active]))
                df = NA
            }
            else
                se = df = 0
            c(est, se, df)
        }
        else c(NA,NA,NA)
    })
    t(result)
}

.clm.vcovHook = function(object, tol = 1e-8, ...) {
    scols = object@misc$scale.idx
    bhat = object@bhat
    active = !is.na(bhat)
    bhat[!active] = 0
    linfct = object@linfct
    rsigma = as.numeric(linfct[, scols, drop = FALSE] %*% object@bhat[scols])
    rsigma = exp(rsigma) * apply(linfct, 1, .is.estble, object@nbasis, tol)
    # I'll do the scaling later
    eta = as.numeric(linfct[, -scols, drop = FALSE] %*% bhat[-scols])
    if (!is.null(object@grid$.offset.))
        eta = eta + object@grid$.offset.
    for (j in scols) linfct[, j] = eta * linfct[, j]
    linfct = (diag(rsigma) %*% linfct) [, active, drop = FALSE]
    linfct %*% tcrossprod(object@V, linfct)
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
