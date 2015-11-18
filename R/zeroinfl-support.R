# Preliminary Support for ZIP models

# zeroinfl objects (pscr package)

recover.data.zeroinfl = function(object, mode = c("response", "count", "zero"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode %in% c("count", "zero"))
        trms = delete.response(terms(object, model = mode))
    else {
        trms = delete.response(terms(object, model = "count"))
        zv = all.vars(delete.response(terms(object, model = "zero")))
        if (length(zv) > 0)
            trms = update(trms, reformulate(zv))
    }
    recover.data(fcall, trms, object$na.action, ...)
}


lsm.basis.zeroinfl = function(object, trms, xlev, grid, 
        mode = c("response", "count", "zero"), ...) 
{
    mode = match.arg(mode)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    if (mode %in% c("count", "zero")) {
        contr = object$contrasts[[mode]]
        X = model.matrix(trms, m, contrasts.arg = contr)
        bhat = coef(object, model = mode)
        V = vcov(object, model = mode)
        if (mode == "count")
            misc = list(tran = "log", inv.lbl = "count", estName = "log(count)")
        else
            misc = list(tran = object$link, inv.lbl = "prob", 
                        estName = paste(object$link, "(prob)", sep = ""))
    }
    else {
        trms1 = delete.response(terms(object, model = "count"))
        contr1 = object$contrasts[["count"]]
        X1 = model.matrix(trms1, m, contrasts.arg = contr1)
        b1 = coef(object, model = "count")
        eta1 = as.numeric(exp(X1 %*% b1))
        
        trms2 = delete.response(terms(object, model = "zero"))
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        eta2 = as.numeric(object$linkinv(X2 %*% b2))
        
        X = diag(eta1) %*% cbind(diag(1 - eta2) %*% X1, diag(-eta2) %*% X2)
        V = X %*% tcrossprod(vcov(object, model = "full"), X)
        bhat = (1 - eta2) * eta1
        X = diag(1, length(bhat))
        
        misc = list(estName = "response")
    }
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) object$df.residual
    dfargs = list()
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}



