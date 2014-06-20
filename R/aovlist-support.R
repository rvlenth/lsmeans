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
    colsums = apply(X[, setdiff(xnms, "(Intercept)"), drop=FALSE], 2, sum)
    if (any(round(colsums,3) != 0))
        warning("Some predictors are correlated with the intercept - results are biased.\n",
                "May help to re-fit with different contrasts, e.g. 'contr.sum'")
    message("NOTE: 'aovlist' results use intra-block analyses and are sometimes flaky")
    #     if(!all(unlist(contr) %in% c("contr.sum","contr.helmert","contr.poly")))
    #         warning("Model uses non-orthogonal contrasts - lsmeans results probably incorrect\n",
    #             "Contrast results probably OK if on one stratum (integer df)")
    
    # initialize arrays
    nonint = setdiff(names(object), "(Intercept)")
    
    k = length(xnms)
    bhat = rep(NA, k) # I'll use NAs to track which slots I've filled
    V = matrix(0, nrow=k, ncol=k)
    names(bhat) = xnms
    dimnames(V) = list(xnms, xnms)
    empty.list = as.list(nonint)
    names(empty.list) = nonint
    Vmats = Vidx = Vdf = empty.list
    wts = matrix(0, nrow = length(nonint), ncol = k)
    dimnames(wts) = list(nonint, xnms)
    # NOTE: At present, I just do intra-block analysis: wts are all 0 and 1
    # Some earlier experimental code is commented-out using '#--'
    #--bmat = wts
    btemp = bhat #++ temp for tracking indexes
    #++Work thru strata in reverse order
    for (nm in rev(nonint)) {
        x = object[[nm]]
        bi = coef(x)
        bi = bi[!is.na(bi)]
        ii = match(names(bi), xnms)
        #--bmat[nm,ii] = bi
        #++ stmts for intra-block analysis
        Vidx[[nm]] = use = setdiff(ii, which(!is.na(bhat))) #++ omit elts already filled
        if(length(use) > 0) {
            ii.left = seq_along(ii)[!is.na(match(ii,use))]
            wts[nm, use] = 1
            bhat[use] = bi[ii.left]
            Vi = vcov(x)[ii.left,ii.left]
            Vi[is.nan(Vi)] = 0 # replace NaNs with 0 - no est of error
            Vmats[[nm]] = Vi
            V[use,use] = Vi
        }
        else {
            Vmats[[nm]] = matrix(0, nrow=0, ncol=0)
        }
        #++ end stmts for intra-block analysis
        Vdf[[nm]] = x$df
        #--piv = x$qr$pivot[1L:x$rank]
        #--wts[nm,ii[piv]] = diag(x$qr$qr)[seq_along(ii)]^2
    }
    #--colsums = apply(wts,2,sum)
    #--colsums[zapsmall(colsums) == 0] = 1
    #--for (i in seq_along(nonint))
    #--     wts[i, ] = wts[i, ] / colsums ## These are efficiencies, not weights
    #--bmat = bmat * wts
    #--bhat = apply(bmat, 2, sum)
    
    x <- object[["(Intercept)"]]
    if (!is.null(x)) {
        # The intercept belongs in the 1st error stratum
        # So have to add a row and column to its covariance matrix
        bhat[1] = x$coefficients[1]
        wts[1,1] = 1
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
    dfargs = list(Vmats=Vmats, Vidx=Vidx, Vdf=unlist(Vdf), wts = wts)
    dffun = function(k, dfargs) {
        v = sapply(seq_along(dfargs$Vdf), function(j) {
            ii = dfargs$Vidx[[j]]
            kk = (k * dfargs$wts[j, ])[ii]            
            sum(kk * dfargs$Vmats[[j]] %*% kk)
        })
        sum(v)^2 / sum(v^2 / dfargs$Vdf) # Good ole Satterthwaite
    }
    nbasis = matrix(NA)  # Consider this further?
    misc = list()
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs, misc = misc)
}