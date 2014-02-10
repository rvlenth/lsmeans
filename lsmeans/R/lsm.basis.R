# create basis for the reference grid

# Generic S3 method - args:
#     object - the model object
#     trms   - terms component of object
#     xlev   - named list of factor levels (but not the coerced ones)
#     grid   - reference grid
# All methods must return a list with these elements:
#     X      - basis for linear fcns over reference grid
#     bhat   - regression coefficients for fixed effects (INCLUDING any NAs)
#     nbasis - matrix whose columns for a basis for non-estimable functions of beta; matrix(NA) if none
#     V      - estimated covariance matrix of bhat
#     dffun  - function(k, dfargs) to find df for k'bhat having std error se
#     dfargs  - additional arguments, if any, for dffun
# Note: if no df exists, set dffun = function(...) NA and dfargs = list()
    
lsm.basis = function(object, trms, xlev, grid)
    UseMethod("lsm.basis")


# for lm objects
lsm.basis.lm <- function(object, trms, xlev, grid) {
    contrasts = attr(model.matrix(object), "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    # coef() works right for lm but coef.aov tosses out NAs
    bhat = as.numeric(object$coefficients) # stretches it out if multivariate
    V = vcov(object)
    
    nbasis = matrix(NA)
    if (sum(is.na(bhat)) > 0) { # get null basis from qr decomp
        # null space of X is same as null space of R in QR decomp
        tR = t(qr.R(object$qr))
        if (ncol(tR) < nrow(tR)) # add columns if not square
            tR = cbind(tR, matrix(0, nrow=nrow(tR), ncol=nrow(tR)-ncol(tR)))
        rank = object$qr$rank
        # last few rows are zero -- add a diagonal
        for (i in (rank+1):nrow(tR)) 
            tR[i,i] = 1
        nbasis = qr.resid(qr(tR[, seq_len(rank)]), tR[, -seq_len(rank)])
        if (!is.matrix(nbasis)) 
            nbasis = matrix(nbasis, ncol=1)
        # permute the rows via pivot
        nbasis[object$qr$pivot, ] = nbasis
    }
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}

# Extension for multivariate case
lsm.basis.mlm <- function(object, trms, xlev, grid) {
    bas <- lsm.basis.lm(object, trms, xlev, grid)
    bhat = coef(object)
    k = ncol(bhat)
    bas$X = kronecker(diag(rep(1,k)), bas$X)
    bas$nbasis = kronecker(rep(1,k), bas$nbasis)
    ylevs = dimnames(bhat)[[2]]
    if (is.null(ylevs)) ylevs = 1:k
    # Quirky, but I use dfargs as repository for ylevs
    bas$dfargs$ylevs = list(rep.meas = ylevs)
    bas
}

lsm.basis.merMod <- function(object, trms, xlev, grid) {
    bhat = fixef(object)
    contrasts = attr(model.matrix(object), "contrasts")
    V = as.matrix(vcov(object))
    dfargs = list()
    if (isLMM(object)) {
        if (require("pbkrtest")) {
            dfargs = list(unadjV = V, adjV = vcovAdj(object, 0))
            V = as.matrix(dfargs$adjV)
            dffun = function(k, dfargs) .KRdf.mer (dfargs$adjV, dfargs$unadjV, k)
        }
        else {
            warning("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
            dffun = function(k, dfargs) NA
        }
    }
    else if (isGLMM(object))
        dffun = function(k, dfargs) NA
    else 
        stop("Can't handle a nonlinear mixed model")
    
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    list(X=X, bhat=bhat, nbasis=matrix(NA), V=V, dffun=dffun, dfargs=dfargs)
}

# For lme4.0, I think
lsm.basis.mer = function(object, trms, xlev, grid)
    lsm.basis.merMod(object,  trms, xlev, grid)

lsm.basis.lme <- function(object, trms, xlev, grid) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)
    V = vcov(object)
    nbasis = matrix(NA)
    dffun = function(...) NA
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=list())
}

lsm.basis.gls <- function(object, trms, xlev, grid) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = coef(object)
    V = vcov(object)
    nbasis = matrix(NA)
    dfargs = list(df = object$dims$N - object$dims$p)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}


# ----------- addition for polr objects in MASS package ------------
recover.data.polr <- function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)))
}

lsm.basis.polr <- function(object, trms, xlev, grid) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    # Strip out the intercept (borrowed code from predict.polr)
    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) 
        X <- X[, -xint, drop = FALSE]
    bhat = c(coef(object), object$zeta)
    V = vcov(object)
    k = length(object$zeta)
    j = matrix(1, nrow=k, ncol=1)
    J = matrix(1, nrow=nrow(X), ncol=1)
    X = cbind(kronecker(j, X), kronecker(diag(1,k), J))
    dfargs = list(ylevs = list(cut = names(object$zeta)))
    nbasis = matrix(NA)
    dffun = function(...) NA
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}


# default method is an error condition...
lsm.basis.default <- function(object, trms, xlev, grid) {
    stop("Can't handle an object of class ", dQuote(class(object)[1]), "\n",
         .show_supported())
}


# get a list of supported objects
# does this by looking in namespace [ns] and methods [meth]
# then strips that off leaving extensions
.show_supported = function(ns = "lsmeans", meth = "lsm.basis") {
    pat = paste(meth, ".", sep="")
    objs = ls(envir = getNamespace(ns), pat = pat)
    clss = gsub(pat, "", objs)
    c("Objects of the following classes are supported:\n",
      paste(dQuote(setdiff(clss, "default")), collapse = ", "))
}
