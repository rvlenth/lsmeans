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

lsm.basis.clm = function (object, trms, xlev, grid, latent = TRUE, ...) {
    # general stuff
    bhat = coef(object)
    V = vcov(object)
    tJac = object$tJac
    dffun = function(...) NA
    link = as.character(object$info$link)
    nbasis = matrix(NA)
    cnm = dimnames(object$tJac)[[1]]
    if (is.null(cnm))
        cnm = paste(seq_len(nrow(tJac)), "|", 1 + seq_len(nrow(tJac)), sep = "")
    misc = list(ylevs = list(cut = cnm), tran = link, inv.lbl = "cumprob",
                offset.mult = -1)
    
    
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
    
    
    ### ----- Alt return when latent == TRUE ----- ###
    if (latent) {
        # Create constant columns for nominal part of linfct
        nomm = apply(bigNom, 2, mean)
        J = matrix(1, nrow = nrow(X))
        bigX = cbind(kronecker(-J, matrix(nomm, nrow = 1)), X)
        # make zero columns for any scale components
        if ((nscols <- length(bhat) - ncol(bigX)) > 0) {
            zer = matrix(0, ncol = nscols)
            bigX = cbind(bigX, kronecker(J, zer))
        }
        dimnames(bigX)[[2]] = names(bhat)
        return(
            list(X = bigX, bhat = bhat, nbasis = nbasis, V = V, 
                 dffun = dffun, dfargs = list(), misc = list())
        )
    }
    
    # (implied else  !latent)

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
        X = cbind(X, S)
        misc$scale.idx = length(object$alpha) + length(object$beta) + seq_len(ncol(S))
        misc$estHook = as.name(".clm.estHook")
    }
    
    
    ### ----- Piece it together ----- ###
    J = matrix(1, nrow=nrow(tJac))
    bigX = cbind(bigNom, -kronecker(J, X))
    dimnames(bigX)[[2]] = names(bhat)
    
    list(X = bigX, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}

# replacement estimation routine for cases with a scale param
.clm.estHook = function(object, do.se = TRUE, tol = 1e-8) {
    scols = object@misc$scale.idx
    S = object@linfct[, scols, drop = FALSE] # This is actually -S
    linfct = object@linfct
    
    active = which(!is.na(object@bhat))
    bhat = object@bhat[active]
    if (is.null(tol)) 
        tol = 1e-8
    
    result = sapply(seq_len(nrow(linfct)), function(i) {
        if (.is.estble(x, object@nbasis, tol)) {
            rsigma = exp(sum(S[i, ] * object@bhat[scols]))
            x = linfct[i, active] * rsigma
            est = sum(bhat[-scols] * x[-scols])
            if (!is.null(object@grid$.offset.))
                est = est + rsigma * object@grid$.offset.[i] # offset is already negated via offset.mult
            if(do.se) {
                x[scols] = est * x[scols] / rsigma
                se = sqrt(.qf.non0(object@V, x))
                df = NA
            }
            else # if these unasked-for results are used, we're bound to get an error!
                se = df = 0
            c(est, se, df)
        }
        else c(NA,NA,NA)
    })
    t(result)
}


lsm.basis.clmm = function (object, trms, xlev, grid, latent = TRUE, ...) {
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
