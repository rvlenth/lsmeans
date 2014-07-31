### Helper functions for lsmeans
### Here we have 'recover.data' and 'lsm.basis' methods
### For models that this package supports.
#----------------------------------------------------------
### Recover data methods will return a data.frame with 
### the original data, and at least these additional attrs:
#   attr(, "terms")      - terms component of object
#   attr(, "responses")  - names of response variable
#   attr(, "predictors") - names of predictors
#-----------------------------------------------------------
# generic version
recover.data <- function(object, ...)
    UseMethod("recover.data")

#----------------------------------------------------------
### lsm.basis methods create a basis for the reference grid
#
# Required args:
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
#     dfargs - additional arguments, if any, for dffun
#     misc   - Extra info ...
#              -- if extra levels need to be added (e.g. mlm, polr),
#                 put them in misc$ylevs
#              -- For transformations or link fcns, use misc$tran
#                 for name (see 'make.link'), and use misc$inv.lbl
#                 for label to use in 'summary' when tran is inverted
#                 (ref.grid looks at lhs of model for tran if none found)
# Note: if no df exists, set dffun = function(...) NA and dfargs = list()
#--------------------------------------------------------------
# generic version
lsm.basis = function(object, trms, xlev, grid)
    UseMethod("lsm.basis")



#--------------------------------------------------------------
### DEFAULT METHODS (we hit these when a model is NOT supported)
# I'll have it return the message if we caught the error in this way
# Then caller can use try() to check for other types of errors,
# and just print this message otherwise 
recover.data.default <- function(object, ...) {
    paste("Can't handle an object of class ", dQuote(class(object)[1]), "\n",
         paste(.show_supported(), collapse=""))
}

lsm.basis.default <- function(object, trms, xlev, grid) {
    stop("Can't handle an object of class", dQuote(class(object)[1]), "\n",
         .show_supported())
}

# Private fcn to get a list of supported objects
# does this by looking in namespace [ns] and methods [meth]
# then strips that off leaving extensions
.show_supported = function(ns = "lsmeans", meth = "lsm.basis") {
    pat = paste(meth, ".", sep="")
    objs = ls(envir = getNamespace(ns), pattern = pat)
    clss = gsub(pat, "", objs)
    c("Objects of the following classes are supported:\n",
      paste(dQuote(setdiff(clss, "default")), collapse = ", "))
}


#--------------------------------------------------------------
### call' objects
# This recover.data method serves as the workhorse for the others
# For model objects, call this with the object's call and its terms component
# Late addition: if data is non-null, use it in place of recovered data
# Later addition: na.action arg req'd - vector of row indexes removed due to NAs
recover.data.call <- function(object, trms, na.action, data, ...) {
    fcall <- object # because I'm easily confused
    vars <- all.vars(trms)
    tbl = data
    if (is.null(tbl)) {
#         m <- match(c("formula", "data", "subset", "weights", 
#                      "na.action", "offset"), names(fcall), 0L)
# I think we don't need to match some of these just to recover the data        
        m <- match(c("formula", "data", "subset"), names(fcall), 0L)
        fcall <- fcall[c(1L, m)]
        fcall$drop.unused.levels <- TRUE
        fcall[[1L]] <- as.name("model.frame")
        fcall$xlev <- NULL # we'll ignore xlev
        fcall$na.action <- na.omit
        # (moved earlier)  vars <- all.vars(trms) # (length will always be >= 2)
        # Put one var on left - keeps out lhs transformations
        if (length(vars) > 1) 
            form <- reformulate(vars[-1], response = vars[1])
        else 
            form <- reformulate(vars)
        fcall$formula <- update(trms, form)
        env <- environment(trms)
        if (is.null(env)) 
            env <- parent.frame()
        tbl <- eval(fcall, env, parent.frame())
        # Drop rows associated with NAs in data
        if (!is.null(na.action))
            tbl = tbl[-(na.action),  , drop=FALSE]
    }
    
    else
        fcall$data = tbl[complete.cases(data), , drop=FALSE]
    
    attr(tbl, "call") = object # the original call
    attr(tbl, "terms") = trms
    attr(tbl, "predictors") = all.vars(delete.response(trms))
    attr(tbl, "responses") = setdiff(vars, attr(tbl, "predictors"))
    tbl
}


#--------------------------------------------------------------
### lm objects (and also aov, rlm, others that inherit) -- but NOT aovList
recover.data.lm <- function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(terms(object)), object$na.action, ...)
}

lsm.basis.lm <- function(object, trms, xlev, grid) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    # coef() works right for lm but coef.aov tosses out NAs
    bhat = as.numeric(object$coefficients) 
    # stretches it out if multivariate - see mlm method
    V = vcov(object)
    
    if (sum(is.na(bhat)) > 0)
        nbasis = nonest.basis(object$qr)
    else
        nbasis = matrix(NA)
    misc = list()
    if (inherits(object, "glm")) {
        fam = object$family
        misc$tran = fam$link
        misc$inv.lbl = "response"
        if (length(grep("binomial", fam$family)) == 1)
            misc$inv.lbl = "prob"
        else if (length(grep("poisson", fam$family)) == 1)
            misc$inv.lbl = "rate"
        dffun = function(k, dfargs) NA
        dfargs = list()
    }
    else {
        dfargs = list(df = object$df.residual)
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
### mlm objects
# (recover.data.lm works just fine)

lsm.basis.mlm <- function(object, trms, xlev, grid) {
    bas <- lsm.basis.lm(object, trms, xlev, grid)
    bhat = coef(object)
    k = ncol(bhat)
    bas$X = kronecker(diag(rep(1,k)), bas$X)
    bas$nbasis = kronecker(rep(1,k), bas$nbasis)
    ylevs = dimnames(bhat)[[2]]
    if (is.null(ylevs)) ylevs = seq_len(k)
    bas$misc$ylevs = list(rep.meas = ylevs)
    bas
}



#--------------------------------------------------------------
### merMod objects (lme4 package)
recover.data.merMod <- function(object, ...) {
    if(!isLMM(object) && !isGLMM(object)) 
        stop("Can't handle a nonlinear mixed model")
    fcall = object@call
    recover.data(fcall, delete.response(terms(object)), 
                 attr(object@frame, "na.action"), ...)
}

lsm.basis.merMod <- function(object, trms, xlev, grid) {
    V = as.matrix(vcov(object))
    dfargs = misc = list()
    if (isLMM(object)) {
        pbdis = .lsm.is.true("disable.pbkrtest")
        if (!pbdis && requireNamespace("pbkrtest")) {
            dfargs = list(unadjV = V, 
                adjV = pbkrtest::vcovAdj.lmerMod(object, 0))
            V = as.matrix(dfargs$adjV)
            dffun = function(k, dfargs) .KRdf.mer (dfargs$adjV, dfargs$unadjV, k)
        }
        else {
            if(!pbdis) message("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
            dffun = function(k, dfargs) NA
        }
    }
    else if (isGLMM(object)) {
        dffun = function(k, dfargs) NA
        fam = family(object)
        misc$tran = fam$link
        misc$inv.lbl = "response"
        if (length(grep("binomial", fam$family)) == 1)
            misc$inv.lbl = "prob"
        else if (length(grep("poisson", fam$family)) == 1)
            misc$inv.lbl = "rate"
    }
    else 
        stop("Can't handle a nonlinear mixed model")
    
    contrasts = attr(object@pp$X, "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)
    
    if (length(bhat) < ncol(X)) {
        # Newer versions of lmer can handle rank deficiency, but we need to do a couple of
        # backflips to put the pieces together right,
        # First, figure out which columns were retained
        kept = match(names(bhat), dimnames(X)[[2]])
        # Now re-do bhat with NAs in the right places
        bhat = NA * X[1, ]
        bhat[kept] = fixef(object)
        # we have to reconstruct the model matrix
        modmat = model.matrix(trms, object@frame, contrasts.arg=contrasts)
        nbasis = nonest.basis(modmat)
    }
    else
        nbasis=matrix(NA)
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
### mer objects (from old lme4 version, now lme4.0 I guess)
recover.data.mer <- recover.data.merMod

lsm.basis.mer <- lsm.basis.merMod
# Seems OK just now (Apr 2014), as vcovAdj methods are identical
# but I haven't tested this



#--------------------------------------------------------------
### lme objects (nlme package)
recover.data.lme <- recover.data.lm

lsm.basis.lme <- function(object, trms, xlev, grid) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = fixef(object)
    V = vcov(object)
    misc = list()
    if (!is.null(object$family)) {
        fam = object$family
        misc$tran = fam$link
        misc$inv.lbl = "response"
        if (length(grep("binomial", fam$family)) == 1)
            misc$inv.lbl = "prob"
        else if (length(grep("poisson", fam$family)) == 1)
            misc$inv.lbl = "rate"
    }
    nbasis = matrix(NA)
    dffun = function(...) NA
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=list(), misc=misc)
}



#--------------------------------------------------------------
### gls objects (nlme package)
recover.data.gls <- function(object, ...) {
    fcall = object$call
    recover.data(fcall, delete.response(getCovariateFormula(object)), 
                 object$na.action, ...)
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
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=list())
}



#--------------------------------------------------------------
### polr objects (MASS package)
recover.data.polr <- recover.data.lm

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
    # Tricky, tricky: need to reverse the sign of the X part
    # because lin. pred is zeta - eta
    X = cbind(kronecker(-j, X), kronecker(diag(1,k), J))
    link = object$method
    if (link == "logistic") link = "logit"
    misc = list(ylevs = list(cut = names(object$zeta)), 
                tran = link, inv.lbl = "cumprob")
    nbasis = matrix(NA)
    dffun = function(...) NA
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=list(), misc=misc)
}




#--------------------------------------------------------------
### survreg objects (survival package)
recover.data.survreg <- recover.data.lm

# Seems to work right in a little testing.
# However, it fails sometimes if I update the model 
# with a subset argument. Workaround: just fitting a new model
lsm.basis.survreg <- function(object, trms, xlev, grid) {
    # Much of this code is adapted from predict.survreg
    bhat = object$coefficients
    k = length(bhat)
    V = vcov(object)[seq_len(k), seq_len(k), drop=FALSE]
    is.fixeds = (k == ncol(object$var))
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)    
    # Hmmm, differs from my lm method using model.matrix(trms, m, contrasts)
    X = model.matrix(object, m)
    nbasis = nonest.basis(model.matrix(object))
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    if (object$dist %in% c("exponential","weibull","loglogistic","loggaussian","lognormal")) 
        misc = list(tran = "log", inv.lbl = "response")
    else 
        misc = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
###  coxph objects (survival package)
recover.data.coxph <- recover.data.survreg

lsm.basis.coxph <- function(object, trms, xlev, grid) {
    object$dist = "doesn't matter"
    result = lsm.basis.survreg(object, trms, xlev, grid)
    result$dfargs$df = NA
    # mimic code for reference = "sample" in predict.coxph
    result$X = result$X - rep(object$means, each = nrow(result$X))
    result$misc$tran = "log"
    result$misc$inv.lbl = "hazard"
    result
}

# Note: Very brief experimentation suggests coxph.penal also works.
# This is an extension of coxph


#--------------------------------------------------------------
###  coxme objects ####
recover.data.coxme <- recover.data.coxph

# I guess this works because it's based on lme code
lsm.basis.coxme <- function(object, trms, xlev, grid) {
    result <- lsm.basis.lme(object, trms, xlev, grid)
    result$misc$tran = "log"
    result$misc$inv.lbl = "hazard"
    result
}


#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------
### Public utility to use for obtaining an orthonormal basis nor nonestimable functions
# Call with its QR decomp (LAPACK=FALSE), if available
nonest.basis <- function(qrX) {
    if (!is.qr(qrX))
        qrX = qr(qrX, LAPACK=FALSE)
    rank = qrX$rank
    tR = t(qr.R(qrX))
    p = nrow(tR)
    if (rank == p)
        return (matrix(NA))
    # null space of X is same as null space of R in QR decomp
    if (ncol(tR) < p) # add columns if not square
        tR = cbind(tR, matrix(0, nrow=p, ncol=p-ncol(tR)))
    # last few rows are zero -- add a diagonal of 1s
    extras = rank + seq_len(p - rank)
    for (j in extras) tR[j,j] = 1
    # nbasis is last p - rank cols of Q in QR decomp of tR
    nbasis = qr.Q(qr(tR))[ , extras, drop = FALSE]
    # permute the rows via pivot
    nbasis[qrX$pivot, ] = nbasis
    nbasis
}
