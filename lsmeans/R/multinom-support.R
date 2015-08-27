### Multinomial modeling

### Example for testing
### From: http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
# library(foreign)
# ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
# library(nnet)
# ml$prog2 <- relevel(ml$prog, ref = "academic")
# test <- multinom(prog2 ~ ses + write, data = ml)
# 
# ### Make order of vcov rows/cols match order of as.numeric(coefs)
# bhat = coef(test)
# ord = as.numeric(matrix(seq_along(bhat), nrow = nrow(bhat), byrow = TRUE))
# V = vcov(test)[ord, ord]

# same as recover.data.lm
recover.data.multinom = function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)), object$na.action, ...)
}

lsm.basis.multinom = function(object, trms, xlev, grid, 
                              mode = c("latent", "prob"), ...) {
# currently mode is ignored - we just support "latent"    
    bhat = t(coef(object))
    V = .my.vcov(object, ...)
    k = ncol(bhat)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    # recenter for latent predictions
    pat = (rbind(0, diag(k + 1, k)) - 1) / (k + 1)
    X = kronecker(pat, X)
    nbasis = estimability::all.estble
    nbasis = kronecker(rep(1,k), nbasis)
    misc = list(tran = "log", inv.lbl = "e^y")
    dfargs = list(df = object$edf)
    dffun = function(k, dfargs) dfargs$df
    ylevs = list(class = object$lev)
    if (is.null(ylevs)) ylevs = list(class = seq_len(k))
    names(ylevs) = as.character(object$call$formula[[2]])
    misc$ylevs = ylevs
    list(X = X, bhat = as.numeric(bhat), nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}
