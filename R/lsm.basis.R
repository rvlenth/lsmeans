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
#     ddfm   - function(k, se, V, misc) to find df for k'bhat having std error se
#     misc   - additional parameters needed by ddfm
# Note: if no df exists, set ddfm = function(...) NA and misc = list()
    
lsm.basis = function(object, trms, xlev, grid, ...)
    UseMethod("lsm.basis")


# for lm objects
lsm.basis.lm <- function(object, trms, xlev, grid) {
    contrasts = attr(model.matrix(object), "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = as.numeric(coef(object)) # stretches it out if multivariate
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
    misc = list(df = nrow(X) - object$rank)
    ddfm = function(k, se, V, misc) misc$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, ddfm=ddfm, misc=misc)
}

# Extension for multivariate case
lsm.basis.mlm <- function(object, trms, xlev, grid) {
    bas <- lsm.basis.lm(object, trms, xlev, grid)
    bhat = coef(object)
    k = ncol(bhat)
    bas$X = kronecker(diag(rep(1,k)), bas$X)
    bas$nbasis = kronecker(rep(1,k), bas$nbasis)
    ylevs = dimnames(bhat)[[2]]
    if (is.null(ynames)) ylevs = 1:k
    bas$misc$ylevs = ylevs
    bas
}

lsm.basis.merMod <- function(object, trms, xlev, grid) {
    bhat = fixef(object)
    contrasts = attr(model.matrix(object), "contrasts")
    misc = list()
    ddfm = function(...) NA
    V = as.matrix(vcov(object))
    if (isLMM(object)) {
        if (require("pbkrtest")) {
            misc = list(unadjV = V)
            V = as.matrix(vcovAdj(object, 0))
            ddfm = function(k, se, V, misc) .KRdf.mer (V, misc$unadjV, k, se*se)
        }
        else {
            warning("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
        }
    }
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    list(X=X, bhat=bhat, nbasis=matrix(NA), V=V, ddfm=ddfm, misc=misc)
}

