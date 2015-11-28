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
    attr(trms, ".Environment")$offset = stats::offset
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
        eta1 = as.numeric(exp(X1 %*% b1 + off1))
        
        trms2 = delete.response(terms(object, model = "zero"))
        off2 = .get.offset(trms2, grid)
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        lp2 = as.numeric(X2 %*% b2)
        eta2 = object$linkinv(lp2 + off2)
        eta2prime = stats::make.link(object$link)$mu.eta(lp2)
        
        delta = diag(eta1) %*% cbind(diag(1 - eta2) %*% X1, diag(-eta2prime) %*% X2)
        V = delta %*% tcrossprod(.pscl.vcov(object, model = "full", ...), delta)
        bhat = (1 - eta2) * eta1
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
    attr(trms, ".Environment")$offset = stats::offset
    recover.data(fcall, trms, object$na.action, ...)
}

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
        eta1 = as.numeric(exp(X1 %*% b1 + off1))
        theta1 = object$theta["count"]
        p1 = switch(object$dist$count,
                poisson = 1 - exp(-eta1),
                negbin = 1 - (theta1 / (eta1 + theta1))^theta1,
                geometric = theta1 / (1 + theta1)  )
        dp1 = switch(object$dist$count,
                poisson = exp(-eta1),
                negbin = (theta1 / (eta1 + theta1))^(1 + theta1),
                geometric = 1 / (1 + theta1)^2  )
        
        trms2 = delete.response(terms(object, model = "zero"))
        off2 = .get.offset(trms2, grid)
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        linv = switch(object$dist$zero, binomial = object$linkinv, exp)
        eta2 = as.numeric(linv(X2 %*% b2 + off2))
        theta2 = object$theta["zero"]
        p2 = switch(object$dist$zero,
                binomial = eta2,
                poisson = 1 - exp(-eta2),
                negbin = 1 - (theta2 / (eta2 + theta2))^theta2,
                geometric = theta2 / (1 + theta2)  )
        dp2 = switch(object$dist$zero,
                binomial = 1,
                poisson = exp(-eta2),
                negbin = (theta2 /(eta2 + theta2))^(1 + theta2),
                geometric = 1 / (1 + theta2)^2)

        if (mode == "mean") {
            bhat = p2 * eta1 / p1
            delta = cbind(diag(bhat*(1 - eta1 * dp1 / p1)) %*% X1,
                          diag(eta1 * dp2 / p1) %*% X2)
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
