# Support for zeroinfl and hurdle models

# zeroinfl objects (pscl package)

recover.data.zeroinfl = function(object, mode = c("mean", "count", "zero"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode %in% c("count", "zero"))
        trms = delete.response(terms(object, model = mode))
    else ### mode = "mean"
        trms = delete.response(object$terms$full)
    # following seems to be needed in order for offsets to be supported
    # attr(trms, ".Environment")$offset = stats::offset
    recover.data(fcall, trms, object$na.action, ...)
}


lsm.basis.zeroinfl = function(object, trms, xlev, grid, 
        mode = c("mean", "count", "zero"), ...) 
{
    mode = match.arg(mode)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    if (mode %in% c("count", "zero")) {
        contr = object$contrasts[[mode]]
        X = model.matrix(trms, m, contrasts.arg = contr)
        bhat = coef(object, model = mode)
        V = .pscl.vcov(object, model = mode, ...)
        if (mode == "count")
            misc = list(tran = "log", inv.lbl = "count")
        else
            misc = list(tran = object$link, inv.lbl = "prob")
    }
    else {   ### mode = "mean"
        trms1 = delete.response(terms(object, model = "count"))
        off1 = .get.offset(trms1, grid)
        contr1 = object$contrasts[["count"]]
        X1 = model.matrix(trms1, m, contrasts.arg = contr1)
        b1 = coef(object, model = "count")
        mu1 = as.numeric(exp(X1 %*% b1 + off1))
        
        trms2 = delete.response(terms(object, model = "zero"))
        off2 = .get.offset(trms2, grid)
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        lp2 = as.numeric(X2 %*% b2) + off2
        mu2 = object$linkinv(lp2)
        mu2prime = stats::make.link(object$link)$mu.eta(lp2)
        
        delta = diag(mu1) %*% cbind(diag(1 - mu2) %*% X1, diag(-mu2prime) %*% X2)
        V = delta %*% tcrossprod(.pscl.vcov(object, model = "full", ...), delta)
        bhat = (1 - mu2) * mu1
        X = diag(1, length(bhat))
        
        misc = list(estName = "response", offset.mult = 0)
    }
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) object$df.residual
    dfargs = list()
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}



#### Support for hurdle models

recover.data.hurdle = function(object, mode = c("mean", "count", "zero", "prob.ratio"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode %in% c("count", "zero"))
        trms = delete.response(terms(object, model = mode))
    else ### mode = "mean" or "prob.ratio"
        trms = delete.response(object$terms$full)
    # needed in order for offsets to be supported
    # attr(trms, ".Environment")$offset = stats::offset
    recover.data(fcall, trms, object$na.action, ...)
}

# see expl notes afterward for notations in some of this
lsm.basis.hurdle = function(object, trms, xlev, grid, 
                              mode = c("mean", "count", "zero", "prob.ratio"), ...) 
{
    mode = match.arg(mode)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    if (mode %in% c("count", "zero")) {
        contr = object$contrasts[[mode]]
        X = model.matrix(trms, m, contrasts.arg = contr)
        bhat = coef(object, model = mode)
        V = .pscl.vcov(object, model = mode, ...)
        misc = switch(object$dist[[mode]],
                        binomial = list(tran = object$link, inv.lbl = "prob"),
                        list(tran = "log", inv.lbl = "count"))
     }
    else {   ### mode %in% c("mean", "prob.ratio")
        trms1 = delete.response(terms(object, model = "count"))
        off1 = .get.offset(trms1, grid)
        contr1 = object$contrasts[["count"]]
        X1 = model.matrix(trms1, m, contrasts.arg = contr1)
        b1 = coef(object, model = "count")
        mu1 = as.numeric(exp(X1 %*% b1 + off1))
        theta1 = object$theta["count"]
        p1 = switch(object$dist$count,
                poisson = 1 - exp(-mu1),
                negbin = 1 - (theta1 / (mu1 + theta1))^theta1,
                geometric = 1 - 1 / (1 + mu1)  )
        dp1 = switch(object$dist$count,
                poisson = mu1 * exp(-mu1),
                negbin = mu1 * (theta1 / (mu1 + theta1))^(1 + theta1),
                geometric = mu1 / (1 + mu1)^2  )
        
        trms2 = delete.response(terms(object, model = "zero"))
        off2 = .get.offset(trms2, grid)
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        lp2 = as.numeric(X2 %*% b2 + off2)
        mu2 = switch(object$dist$zero, 
                     binomial = object$linkinv(lp2),
                     exp(lp2)  )
        theta2 = object$theta["zero"]
        p2 = switch(object$dist$zero,
                binomial = mu2,
                poisson = 1 - exp(-mu2),
                negbin = 1 - (theta2 / (mu2 + theta2))^theta2,
                geometric = 1 - 1 / (1 + mu2)  )
        dp2 = switch(object$dist$zero,
                binomial = stats::make.link(object$link)$mu.eta(lp2),
                poisson = mu2 * exp(-mu2),
                negbin = mu2 * (theta2 /(mu2 + theta2))^(1 + theta2),
                geometric = mu2 / (1 + mu2)^2  )

        if (mode == "mean") {
            bhat = p2 * mu1 / p1
            delta = cbind(diag(bhat*(1 - mu1 * dp1 / p1)) %*% X1,
                          diag(mu1 * dp2 / p1) %*% X2)
        }
        else {  ## mode == "prob.ratio"
            bhat = p2 / p1
            delta = cbind(diag(-p2 * dp1 / p1^2) %*% X1,
                          diag(dp2 / p1) %*% X2)
        }
        V = delta %*% tcrossprod(.pscl.vcov(object, model = "full", ...), delta)
        X = diag(1, length(bhat))
        
        misc = list(estName = mode, offset.mult = 0)
    }
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) object$df.residual
    dfargs = list()
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}


# special version of .my.vcov that accepts (and requires!) model argument
.pscl.vcov = function(object, model, vcov. = stats::vcov, ...) {
    if (is.function(vcov.))
        vcov. = vcov.(object, model = model)
    else if (!is.matrix(vcov.))
        stop("vcov. must be a function or a square matrix")
    vcov.
}

# Explanatory notes for hurdle models
# -----------------------------------
#     We have a linear predictor eta = X%*%beta + offset
#     mu = h(eta) where h is inverse link (usually exp but not always)
#     Define p = P(Y > 0 | mu). This comes out to...
#         binomial: mu
#         poisson: 1 - exp(-mu)
#         negbin: 1 - (theta/(mu+theta))^theta
#         geometric: 1 - 1/(mu+1)
#     Define dp = dp/d(eta). Note - when h(mu)=exp(mu) we have dp = mu*dp/d(mu)
#         binomial: h'(eta)
#         poisson: mu*exp(-mu)
#         negbin: mu*(theta/(mu+theta))^(theta+1)
#         geometric: mu/(mu+1)^2
#     
#     This gives us what we need to find the estimates and apply the delta method
#     In the code we index these notations with 1 (count model) and 2 (zero model)
#     And we treat theta1 and theta2 as constants
#
#!!! In theory, above seems correct, and estimates match those from predict.hurdle.
#!!! But SEs don't seem right. 
#!!! They do seem right though if I omit the factor of mu in dp
#!!! when link is log
