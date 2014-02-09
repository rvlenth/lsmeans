# lsmeans support for 'survival' package

# (identical to .lm method)
recover.data.survreg <- recover.data.lm


### Seems to work right in a little testing.
### However, it fails if I update the model to fewer blocks
lsm.basis.survreg <- function(object, trms, xlev, grid) {
    # Much of this code is adapted from predict.survreg
    bhat = object$coefficients
    k = length(bhat)
    V = vcov(object)[1:k, 1:k, drop=FALSE]
    is.fixeds = (k == ncol(object$var))
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)    
    # Hmmm, differs from my lm method using model.matrix(trms, m, contrasts)
    X = model.matrix(object, m)
#    strata <- rep(1L, nrow(X))
#    nstrata <- 1
#    offset <- 0
#    scale <- object$scale[strata]
#     if (is.character(object$dist)) 
#         dd <- survreg.distributions[[object$dist]]
#     else dd <- object$dist
#     if (is.null(dd$itrans)) {
#         itrans <- function(x) x
#         dtrans <- function(x) 1
#     }
#     else {
#         itrans <- dd$itrans
#         dtrans <- dd$dtrans
#     }
#     if (!is.null(dd$dist)) 
#         dd <- survreg.distributions[[dd$dist]]
#     pred <- drop(X %*% bhat) + offset
    
    nbasis = matrix(NA)
    if (sum(is.na(bhat)) > 0) { # get null basis from qr decomp
        warn("No provision yet made for rank-deficient survreg models.\n",
             "LS means results may be misleading")
    }
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs)
}



### coxph objects
recover.data.coxph <- recover.data.survreg

lsm.basis.coxph <- function(object, trms, xlev, grid) {
    result = lsm.basis.survreg(object, trms, xlev, grid)
    # asymptotic
    result$dfargs$df = NA
    # mimic reference = "sample" in predict.coxph
    result$X = result$X - rep(object$means, each = nrow(result$X))
    result
}


