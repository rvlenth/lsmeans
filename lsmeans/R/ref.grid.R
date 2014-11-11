# Reference grid code


# Change to cov.reduce specification: can be...
#     a function: is applied to all covariates
#     named list of functions: applied to those covariates (else mean is used)
#     TRUE - same as mean
#     FALSE - same as function(x) sort(unique(x))

ref.grid <- function(object, at, cov.reduce = mean, mult.name, mult.levs, 
                     options = lsm.options()$ref.grid, data, ...) {
    # recover the data
    if (missing(data)) {
        data = try(recover.data (object, data = NULL))
        if (inherits(data, "try-error"))
            stop("Possible remedy: Supply the data used in the 'data' argument")
    }
    else # attach needed attributes to given data
        data = recover.data(object, data = data)
    
    if(is.character(data)) # 'data' is in fact an error message
        stop(data)
        
    
    trms = attr(data, "terms")
    
    # find out if any variables are coerced to factors
    ### OLD VERSION: anm = all.names(attr(data, "terms"))    
    ###              coerced = anm[1 + grep("factor|ordered", anm)]
    coerced = .find.coerced(trms, data)
    
    # convenience function
    sort.unique = function(x) sort(unique(x))
    
    # Ensure cov.reduce is a function or list thereof
    dep.x = list() # list of formulas to fit later
    fix.cr = function(cvr) {
        # cvr is TRUE or FALSE
        if(is.logical(cvr)) 
            if(cvr[1]) cvr = mean
        else              cvr = sort.unique
        else if (inherits(cvr, "formula")) {
            if (length(cvr) < 3)
                stop("Formulas in 'cov.reduce' must be two-sided")
            lhs = all.vars(cvr)[1]
            dep.x[[lhs]] <<- cvr
            cvr = mean 
        }
        else if (!inherits(cvr, c("function","list")))
            stop("Invalid 'cov.reduce' argument")
        cvr
    }
    
    # IMPORTANT: following stmts may also affect x.dep
    if (is.list(cov.reduce))
        cov.reduce = lapply(cov.reduce, fix.cr)
    else
        cov.reduce = fix.cr(cov.reduce)
    
    # zap any formulas that are also in 'at'
    if (!missing(at))
        for (xnm in names(at)) dep.x[[xnm]] = NULL
    
    
    # local cov.reduce function that works with function or named list
    cr = function(x, nm) {
        if (is.function(cov.reduce))
            cov.reduce(x)
        else if (!is.null(cov.reduce[[nm]]))
            cov.reduce[[nm]](x)
        else
            mean(x)
    }
    
    # initialize empty lists
    ref.levels = matlevs = xlev = list()
    
    for (nm in attr(data, "responses")) {
        y = data[[nm]]
        if (is.matrix(y))
            matlevs[[nm]] = apply(y, 2, mean)
        else
            ref.levels[[nm]] = mean(y)
    }
    
    for (nm in attr(data, "predictors")) {
        x = data[[nm]]
        
    # Save the original levels of factors, no matter what
    if (is.factor(x))
        xlev[[nm]] = levels(factor(x))
    # (applying factor drops any unused levels)
    
    # Now go thru and find reference levels...
        # mentioned in 'at' list but not coerced
        if (!(nm %in% coerced) && !missing(at) && !is.null(at[[nm]]))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x))
            ref.levels[[nm]] = levels(factor(x))
        # matrices
        else if (is.matrix(x)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matlevs[[nm]] = apply(x, 2, cr, nm)
            # if cov.reduce returns a vector, average its columns
            if (is.matrix(matlevs[[nm]]))
                matlevs[[nm]] = apply(matlevs[[nm]], 2, mean)
        }
        # covariate coerced, or not mentioned in 'at'
        else {
            # single numeric pred but coerced to a factor - use unique values
            # even if in 'at' list. We'll fix this up later
            if (nm %in% coerced)            
                ref.levels[[nm]] = sort.unique(x)
            
            # Ordinary covariates - summarize
            else 
                ref.levels[[nm]] = cr(x, nm)
        }
    }
    
    # Now create the reference grid
    grid = do.call(expand.grid, ref.levels)
    
    # add any matrices
    for (nm in names(matlevs))
        grid[[nm]] = matrix(rep(matlevs[[nm]], each=nrow(grid)), nrow=nrow(grid))

    # resolve any covariate formulas
    for (xnm in names(dep.x)) {
        if (!all(all.vars(dep.x[[xnm]]) %in% names(grid)))
            stop("Formulas in 'cov.reduce' must predict covariates actually in the model")
        xmod = lm(dep.x[[xnm]], data = data)
        grid[[xnm]] = predict(xmod, newdata = grid)
        ref.levels[[xnm]] = NULL
    }
    
    basis = lsm.basis(object, trms, xlev, grid, ...)
    
    misc = basis$misc
    
    form = attr(data, "call")$formula
    if (is.null(misc$tran) && (length(form) > 2)) { # No link fcn, but response may be transformed
        lhs = form[[2]]
        tran = setdiff(all.names(lhs), c(all.vars(lhs), "~", "cbind"))
        if(length(tran) == 1)
            misc$tran = tran
    }
    
    # Take care of multivariate response
    multresp = character(0) ### ??? was list()
    ylevs = misc$ylevs
    if(!is.null(ylevs)) { # have a multivariate situation
       if (missing(mult.levs)) {
            if (missing(mult.name))
                mult.name = names(ylevs)[1]
            ref.levels[[mult.name]] = ylevs[[1]]
            multresp = mult.name
            MF = data.frame(ylevs)
            names(MF) = mult.name
        }
        else {
            k = prod(sapply(mult.levs, length))
            if (k != length(ylevs[[1]])) 
                stop("supplied 'mult.levs' is of different length than that of multivariate response")
            for (nm in names(mult.levs))
                ref.levels[[nm]] = mult.levs[[nm]]
            multresp = names(mult.levs)
            MF = do.call("expand.grid", mult.levs)
        }
        ###grid = do.call("expand.grid", ref.levels)
        grid = merge(grid, MF)
        # add any matrices
        for (nm in names(matlevs))
            grid[[nm]] = matrix(rep(matlevs[[nm]], each=nrow(grid)), nrow=nrow(grid))
    }

# Here's a complication. If a numeric predictor was coerced to a factor, we had to
# include all its levels in the reference grid, even if altered in 'at'
# Moreover, whatever levels are in 'at' must be a subset of the unique values
# So we now need to subset the rows of the grid and linfct based on 'at'

    problems = if (!missing(at)) 
        intersect(c(multresp, coerced), names(at)) 
    else character(0)
    if (length(problems > 0)) {
        incl.flags = rep(TRUE, nrow(grid))
        for (nm in problems) {
            if (is.numeric(ref.levels[[nm]])) {
                at[[nm]] = round(at[[nm]], 3)
                ref.levels[[nm]] = round(ref.levels[[nm]], 3)
            }
            # get only "legal" levels
            at[[nm]] = at[[nm]][at[[nm]] %in% ref.levels[[nm]]]
            # Now which of those are left out?
            excl = setdiff(ref.levels[[nm]], at[[nm]])
            for (x in excl)
                incl.flags[grid[[nm]] == x] = FALSE
            ref.levels[[nm]] = at[[nm]]
        }
        if (!any(incl.flags))
            stop("Reference grid is empty due to mismatched levels in 'at'")
        grid = grid[incl.flags, , drop=FALSE]
        basis$X = basis$X[incl.flags, , drop=FALSE]
    }

    # Any offsets??? (misc$offset.mult might specify removing or reversing the offset)
    if(!is.null(attr(trms,"offset"))) {
        om = 1
        if (!is.null(misc$offset.mult))
            om = misc$offset.mult
        if (any(om != 0))
            grid[[".offset."]] = om * .get.offset(trms, grid)
    }

    ### --- Determine frequencies --- (added ver.2.11)
    nms = union(names(xlev), coerced) # only factors, no covariates or mult.resp
    # originally, I used 'plyr::count', but there are probs when there is a 'freq' variable
    id = plyr::id(data[, nms, drop = FALSE], drop = TRUE)
    uid = !duplicated(id)
    key = do.call(paste, data[uid, nms, drop = FALSE])
    key = key[order(id[uid])]
    frq = tabulate(id, attr(id, "n"))
    tgt = do.call(paste, grid[, nms, drop = FALSE])
    freq = rep(0, nrow(grid))
    for (i in seq_along(key))
        freq[tgt == key[i]] = frq[i] ###ftbl[i, "freq"]
    grid[[".freq."]] = freq

    misc$ylevs = NULL # No longer needed
    misc$estName = "prediction"
    misc$infer = c(FALSE,FALSE)
    misc$level = .95
    misc$adjust = "none"
    misc$famSize = nrow(grid)
    misc$avgd.over = character(0)
    
    result = new ("ref.grid",
         model.info = list(call = attr(data,"call"), terms = trms, xlev = xlev),
         roles = list(predictors = attr(data, "predictors"), 
                      responses = attr(data, "responses"), 
                      multresp = multresp),
         grid = grid, levels = ref.levels, matlevs = matlevs,
         linfct = basis$X, bhat = basis$bhat, nbasis = basis$nbasis, V = basis$V,
         dffun = basis$dffun, dfargs = basis$dfargs, misc = misc)

    if(!is.null(options)) {
        options$object = result
        result = do.call("update.ref.grid", options)
    }

    if(!is.null(hook <- misc$postGridHook)) {
        if (is.character(hook))
            hook = get(hook)
        result@misc$postGridHook = NULL
        hook(result)
    }
    else
        result
}


#### End of ref.grid ------------------------------------------




# This function figures out which covariates in a model 
# have been coerced to factors. Does NOT rely on the names of
# functions like 'factor' or 'interaction' as we use actual results
.find.coerced = function(trms, data) {
    isfac = sapply(data, function(x) inherits(x, "factor"))
    
    # Character vectors of factors and covariates in the data...
    facs.d = names(data)[isfac]
    covs.d = names(data)[!isfac]
    
    lbls = attr(trms, "term.labels")
    M = model.frame(trms, data)
    isfac = sapply(M, function(x) inherits(x, "factor"))
    
    # Character vector of terms in the model frame that are factors ...
    facs.m = names(M)[isfac]
    
    # Exclude the terms that are already factors
    # What's left will be things like "factor(dose)", "interact(dose,treat)", etc
    cterms = setdiff(facs.m, facs.d)
    
    if(length(cterms) == 0) 
        return(cterms)
    # (else) Strip off the function calls
    cvars = lapply(cterms, function(x) all.vars(reformulate(x)))
    
    # Exclude any variables that are already factors
    intersect(unique(unlist(cvars)), covs.d)
}

# calculate the offset for the given grid
.get.offset = function(terms, grid) {
    off.idx = attr(terms, "offset")
    offset = rep(0, nrow(grid))
    tvars = attr(terms, "variables")
    for (i in off.idx)
        offset = offset + eval(tvars[[i+1]], grid)
    offset
}


# Computes the quadratic form y'Xy after subsetting for the nonzero elements of y
.qf.non0 = function(X, y) {
    ii = (zapsmall(y) != 0)
    if (any(ii))
        sum(y[ii] * (X[ii, ii, drop = FALSE] %*% y[ii]))
    else 0
}

# utility to check estimability of x'beta, given nonest.basis
is.estble = function(x, nbasis, tol=1e-8) {
    if (is.matrix(x))
        return(apply(x, 1, is.estble, nbasis, tol))
    if(is.na(nbasis[1]))
        TRUE
    else {
        chk = as.numeric(crossprod(nbasis, x))
        ssqx = sum(x*x) # BEFORE subsetting x
        # If x really small, don't scale chk'chk
        if (ssqx < tol) ssqx = 1
        sum(chk*chk) < tol * ssqx
    }
}

# utility fcn to get est's, std errors, and df
# new arg: do.se -- if FALSE, just do the estimates and return 0 for se and df
# returns a data.frame with an add'l "link" attribute if misc$tran is non-null
# .est.se.df = function(linfct, bhat, nbasis, V, dffun, dfargs, misc, do.se=TRUE, 
#                       tol=getOption("lsmeans")$estble.tol) {
# 2.13: Revised to call w/ just object instead of all those args (except linfct)
# Also moved offest comps to here, and provided for misc$estHook
.est.se.df = function(object, do.se=TRUE, tol=lsm.options()$estble.tol) {
    if (is.null(tol)) 
        tol = 1e-8
    misc = object@misc
    if (!is.null(hook <- misc$estHook)) {
        if (is.character(hook)) hook = get(hook)
        result = hook(object, do.se=do.se, tol=tol)
    }
    else {
        active = which(!is.na(object@bhat))
        bhat = object@bhat[active]
        result = t(apply(object@linfct, 1, function(x) {
            if (is.estble(x, object@nbasis, tol)) {
                x = x[active]
                est = sum(bhat * x)
                if(do.se) {
                    se = sqrt(.qf.non0(object@V, x))
                    df = object@dffun(x, object@dfargs)
                }
                else # if these unasked-for results are used, we're bound to get an error!
                    se = df = 0
                c(est, se, df)
            }
            else c(NA,NA,NA)
        }))
        if (!is.null(object@grid$.offset.))
            result[, 1] = result[, 1] + object@grid$.offset.
    }
    result = as.data.frame(result)
    names(result) = c(misc$estName, "SE", "df")
    
    if (!is.null(misc$tran) && (misc$tran != "none")) {
        if(is.character(misc$tran)) {
            link = try(make.link(misc$tran), silent=TRUE)
            if (!inherits(link, "try-error"))
                attr(result, "link") = link
        }
        else if (is.list(misc$tran))
            attr(result, "link") = misc$tran
    }
    result
}

### =========== Methods for ref.grid class =============================

str.ref.grid <- function(object, ...) {
    showlevs = function(x) { # internal convenience function
        if (is.null(x)) cat("(predicted by other variables)")
        else cat(paste(format(x, digits = 5, justify = "none"), collapse=", "))
    }
    #cat("responses: ")
    #showlevs(object@roles$responses)
    levs = object@levels
    cat(paste("'", class(object)[1], "' object with variables:\n", sep=""))
    for (nm in union(object@roles$predictors, union(object@roles$multresp, object@roles$responses))) {
        cat(paste("    ", nm, " = ", sep = ""))
        if (nm %in% names(object@matlevs)) {
            if (nm %in% object@roles$responses)
                cat("multivariate response with means: ")
            else
                cat("matrix with column means: ")
            cat("\n        ")
            showlevs(object@matlevs[[nm]])
        }
        else if (nm %in% object@roles$multresp) {
            cat("multivariate response levels: ")
            showlevs(levs[[nm]])
        }
        else if (nm %in% object@roles$responses) {
            cat("response variable with mean ")
            showlevs(levs[[nm]])
        }
        else
            showlevs(levs[[nm]])
        cat("\n")
    }
    if(!is.null(tran <- object@misc$tran)) {
        if (is.list(tran)) tran = "custom - see slot(, \"misc\")$tran"
        cat(paste("Transformation:", dQuote(tran), "\n"))
    }
}


# utility to compute an adjusted p value
# tail is -1, 0, 1 for left, two-sided, or right
.adj.p.value = function(t, df, adjust, fam.info, tail) {
# do a pmatch of the adjust method, case insensitive
    adj.meths = c("sidak", "tukey", "scheffe", p.adjust.methods)
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k))
        stop("Adjust method '", adjust, "' is not recognized or not valid")
    adjust = adj.meths[k]
    if ((tail != 0) && (k %in% 2:3)) # One-sided tests, change Tukey/Scheffe to Bonferroni
        adjust = "bonferroni"
    
    # pseudo-asymptotic results when df is NA
    df[is.na(df)] = 10000
    
    fam.size = fam.info[1]
    n.contr = fam.info[2] ## n.contr = sum(!is.na(t))
    abst = abs(t)
    if (tail == 0)
        unadj.p = 2*pt(abst, df, lower.tail=FALSE)
    else
        unadj.p = pt(t, df, lower.tail = (tail<0))
    if (adjust %in% p.adjust.methods) {
        if (n.contr == length(unadj.p))
            pval = p.adjust(unadj.p, adjust, n = n.contr)
        else
            pval = as.numeric(apply(matrix(unadj.p, nrow=n.contr), 2, 
                function(pp) p.adjust(pp, adjust, n=sum(!is.na(pp)))))
    }
    else pval = switch(adjust,
        sidak = 1 - (1 - unadj.p)^n.contr,
        tukey = ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),
        scheffe = pf(t^2/(fam.size-1), fam.size-1, df, lower.tail=FALSE),
    )
    chk.adj = match(adjust, c("none", "tukey", "scheffe"), nomatch = 99)
    do.msg = (chk.adj > 1) && (n.contr > 1) && 
             !((fam.size == 2) && (chk.adj < 10)) 
    if (do.msg) {
        xtra = if(chk.adj < 10) paste("a family of", fam.size, "means")
               else             paste(n.contr, "tests")
        mesg = paste("P value adjustment:", adjust, "method for", xtra)
    }
    else mesg = NULL
    list(pval=pval, mesg=mesg, adjust=adjust)
}

# Code needed for an adjusted critical value
# returns a list similar to .adj.p.value
.adj.critval = function(level, df, adjust, fam.info) {
    mesg = NULL
    adj.meths = c("tukey", "sidak", "bonferroni", "none")
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k)) {
        k = which(adj.meths == "none")
        mesg = "Confidence levels are NOT adjusted for multiplicity"
    }
    adjust = adj.meths[k]
    
    # pseudo-asymptotic results when df is NA
    df[is.na(df)] = 10000
    
    fam.size = fam.info[1]
    n.contr = fam.info[2]
    
    chk.adj = match(adjust, c("none", "tukey"), nomatch = 99)
    do.msg = (chk.adj > 1) && (n.contr > 1) && 
        !((fam.size == 2) && (chk.adj < 10)) 
    if (do.msg) {
        xtra = if(chk.adj < 10) paste("a family of", fam.size, "means")
        else             paste(n.contr, "tests")
        mesg = paste("Confidence-level adjustment:", adjust, "method for", xtra)
    }
    
    cv = switch(adjust,
        none = -qt((1-level)/2, df),
        sidak = -qt((1 - level^(1/n.contr))/2, df),
        bonferroni = -qt((1-level)/n.contr/2, df),
        tukey = qtukey(level, fam.size, df) / sqrt(2)
    )
    list(cv = cv, mesg = mesg, adjust = adjust)
}

### Support for different prediction types ###

# Valid values for type arg or predict.type option
.valid.types = c("link","response","lp","linear")

# get "predict.type" option from misc, and make sure it's legal
.get.predict.type = function(misc) {
    type = misc$predict.type
    if (is.null(type))
        .valid.types[1]
    else
        .validate.type(type)
}

# check a "type" arg to make it legal
.validate.type = function (type) {
    .valid.types[pmatch(type, .valid.types, 1)]
}

# S3 predict method
predict.ref.grid <- function(object, type, ...) {
    # update with any "summary" options
    opt = lsm.options()$summary
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.ref.grid", opt)
    }
    
    if (missing(type))
        type = .get.predict.type(object@misc)
    else
        type = .validate.type(type)
    
    pred = .est.se.df(object, do.se=FALSE)
    result = pred[[1]]
# MOVED TO .EST.SE.DF    
#     if (".offset." %in% names(object@grid))
#         result = result + object@grid[[".offset."]]
    if (type == "response") {
        link = attr(pred, "link")
        if (!is.null(link))
            result = link$linkinv(result)
    }
    result
}

# S3 summary method
summary.ref.grid <- function(object, infer, level, adjust, by, type, df, 
                             null = 0, delta = 0, side = 0, ...) {
    # update with any "summary" options
    opt = lsm.options()$summary
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.ref.grid", opt)
    }
    
    if(missing(df)) df = object@misc$df
    if(!is.null(df))
        object@dffun = function(k, dfargs) df
    
    # reconcile all the different ways we could specify the alternative
    # ... and map each to one of the first 3 subscripts
    side.opts = c("left","both","right","two-sided","noninferiority","nonsuperiority","equivalence","superiority","inferiority","0","2","-1","1","+1","<",">","!=","=")
    side.map =  c( 1,     2,     3,      2,          3,               1,               2,            3,            1,            2,  2,   1,  3,   3,  1,  3,  2,   2)
    side = side.map[pmatch(side, side.opts, 2)[1]] - 2
    delta = abs(delta)
    
    result = .est.se.df(object)
    
    lblnms = setdiff(names(object@grid), 
                     c(object@roles$responses, ".offset.", ".freq."))
    lbls = object@grid[lblnms]

    ### implement my 'variable defaults' scheme    
    if(missing(infer)) infer = object@misc$infer
    if(missing(level)) level = object@misc$level
    if(missing(adjust)) adjust = object@misc$adjust
    if(missing(by)) by = object@misc$by.vars
    
    if (missing(type))
        type = .get.predict.type(object@misc)
    else
        type = .validate.type(type)

    zFlag = (all(is.na(result$df)))
    inv = (type == "response") # flag to inverse-transform
    
    if ((length(infer) == 0) || !is.logical(infer)) 
        infer = c(FALSE, FALSE)
    if(length(infer == 1)) 
        infer = c(infer,infer)

    if(inv && !is.null(object@misc$tran)) {
        link = attr(result, "link")
        if (!is.null(object@misc$inv.lbl))
            names(result)[1] = object@misc$inv.lbl
        else
            names(result)[1] = "lsresponse"
    }
    else
        link = NULL
    
    attr(result, "link") = NULL
    estName = names(result)[1]

    mesg = object@misc$initMesg
    
    by.size = nrow(object@grid)
    if (!is.null(by))
        for (nm in by)
            by.size = by.size / length(unique(object@levels[[nm]]))
    fam.info = c(object@misc$famSize, by.size)
    cnm = NULL
    
    if(infer[1]) { # add CIs
        quant = 1 - (1 - level)/2
        ###cv = if(zFlag) qnorm(quant) else qt(quant, result$df)
        acv = .adj.critval(level, result$df, adjust, fam.info)
        adjust = acv$adjust
        cv = acv$cv
        cnm = if (zFlag) c("asymp.LCL", "asymp.UCL") else c("lower.CL","upper.CL")
        result[[cnm[1]]] = result[[1]] - cv*result$SE
        result[[cnm[2]]] = result[[1]] + cv*result$SE
        if (!is.null(link)) {
            result[[cnm[1]]] = link$linkinv(result[[cnm[1]]])
            result[[cnm[2]]] = link$linkinv(result[[cnm[2]]])
        }
        mesg = c(mesg, paste("Confidence level used:", level), acv$mesg)
    }
    if(infer[2]) { # add tests
        if (!all(null == 0)) {
            result[["null"]] = null
            if (!is.null(link))
                result[["null"]] = link$linkinv(result[["null"]])
        }
        tnm = ifelse (zFlag, "z.ratio", "t.ratio")
        tail = ifelse(side == 0, -sign(abs(delta)), side)
        if (side == 0) {
            if (delta == 0) # two-sided sig test
                t.ratio = result[[tnm]] = (result[[1]] - null) / result$SE
            else
                t.ratio = result[[tnm]] = (abs(result[[1]] - null) - delta) / result$SE
        }
        else {
            t.ratio = result[[tnm]] = (result[[1]] - null + side * delta) / result$SE            
        }
        apv = .adj.p.value(t.ratio, result$df, adjust, fam.info, tail)
        adjust = apv$adjust   # in case it was abbreviated
        result$p.value = apv$pval
        mesg = c(mesg, apv$mesg)
        if (delta > 0)
            mesg = c(mesg, paste("Statistics are tests of", c("nonsuperiority","equivalence","noninferiority")[side+2],
                                 "with a threshold of", delta))
        if(tail != 0) 
            mesg = c(mesg, paste("P values are ", ifelse(tail<0,"left-","right-"),"tailed", sep=""))
        if (!is.null(link)) 
            mesg = c(mesg, "Tests are performed on the linear-predictor scale")
    }
    if (!is.null(link)) {
        result[["SE"]] = link$mu.eta(result[[1]]) * result[["SE"]]
        result[[1]] = link$linkinv(result[[1]])
    }
    
    if (length(object@misc$avgd.over) > 0)
        mesg = c(paste("Results are averaged over the levels of:",
                 paste(object@misc$avgd.over, collapse = ", ")), mesg)

    summ = cbind(lbls, result)
    attr(summ, "estName") = estName
    attr(summ, "clNames") = cnm  # will be NULL if infer[1] is FALSE
    attr(summ, "pri.vars") = setdiff(union(object@misc$pri.vars, object@misc$by.vars), by)
    attr(summ, "by.vars") = by
    attr(summ, "mesg") = unique(mesg)
    class(summ) = c("summary.ref.grid", "data.frame")
    summ
}


# left-or right-justify column labels for m depending on "l" or "R" in just
.just.labs = function(m, just) {
    nm = dimnames(m)[[2]]
    for (j in seq_len(length(nm))) {
        if(just[nm[j]] == "L") 
            nm[j] = format(nm[j], width = nchar(m[1,j]), just="left")
    }
    dimnames(m) = list(rep("", nrow(m)), nm)
    m
}

# Format a data.frame produced by summary.ref.grid
print.summary.ref.grid = function(x, ..., digits=NULL, quote=FALSE, right=TRUE) {
    x.save = x
    if (!is.null(x$df)) x$df = round(x$df, 2)
    if (!is.null(x$t.ratio)) x$t.ratio = round(x$t.ratio, 3)
    if (!is.null(x$p.value)) {
        fp = x$p.value = format(round(x$p.value,4), nsmall=4, sci=FALSE)
        x$p.value[fp=="0.0000"] = "<.0001"
    }
    just = sapply(x.save, function(col) if(is.numeric(col)) "R" else "L")
    xc = as.matrix(format.data.frame(x, digits=digits, na.encode=FALSE))
    m = apply(rbind(just, names(x), xc), 2, function(x) {
        w = max(sapply(x, nchar))
        if (x[1] == "R") format(x[-seq_len(2)], width = w, justify="right")
        else format(x[-seq_len(2)], width = w, justify="left")
    })
    if(!is.matrix(m)) m = t(as.matrix(m))
    by.vars = attr(x, "by.vars")
    if (is.null(by.vars)) {
        m = .just.labs(m, just)
        print(m, quote=FALSE, right=TRUE)
        cat("\n")
    }
    else { # separate listing for each by variable
        m = .just.labs(m[, setdiff(names(x), by.vars)], just)
        pargs = as.list(x[,by.vars, drop=FALSE])
        pargs$sep = ", "
        lbls = do.call(paste, pargs)
        for (lb in unique(lbls)) {
            rows = which(lbls==lb)
            levs = paste(by.vars, "=", xc[rows[1], by.vars])
            cat(paste(paste(levs, collapse=", ")), ":\n", sep="")
            print(m[rows, , drop=FALSE], ..., quote=quote, right=right)
            cat("\n")
        }
    }
    
    msg = unique(attr(x, "mesg"))
    if (!is.null(msg))
        for (j in seq_len(length(msg))) cat(paste(msg[j], "\n"))
    
    invisible(x.save)
}


print.ref.grid = function(x,...)
    print(summary.ref.grid(x, ...))



# vcov method
vcov.ref.grid = function(object, ...) {
    tol = lsm.options()$estble.tol
    if(is.null(tol)) 
        tol = 1e-8
    if (!is.null(hook <- object@misc$vcovHook)) {
        if (is.character(hook)) 
            hook = get(hook)
        hook(object, tol = tol, ...)
    }
    else {
        X = object@linfct
        estble = is.estble(X, object@nbasis, tol) ###apply(X, 1, .is.estble, object@nbasis, tol)
        X[!estble, ] = NA
        X = X[, !is.na(object@bhat)]
        X %*% tcrossprod(object@V, X)
    }
}


# Method to alter contents of misc slot
update.ref.grid = function(object, ..., silent = FALSE) {
    args = list(...)
    valid.choices = c("adjust","alpha","avgd.over","by.vars","df",
        "initMesg","estName","famSize","infer","inv.lbl",
        "level","methdesc","predict.type","pri.vars","tran")
    misc = object@misc
    for (nm in names(args)) {
        fullname = try(match.arg(nm, valid.choices), silent=TRUE)
        if(inherits(fullname, "try-error")) {
            if (!silent)
                message("Argument ", sQuote(nm), " was ignored. Valid choices are:\n",
                    paste(valid.choices, collapse=", "))
        }
        else {
            if (fullname == "by.vars") {
                allvars = union(misc$pri.vars, misc$by.vars)
                misc$pri.vars = setdiff(allvars, args[[nm]])
            }
            if (fullname == "pri.vars") {
                allvars = union(misc$pri.vars, misc$by.vars)
                misc$by.vars = setdiff(allvars, args[[nm]])
            }
            misc[[fullname]] = args[[nm]]
        }
    }
    object@misc = misc
    object
}

### set or change lsmeans options
lsm.options = function(...) {
    opts = getOption("lsmeans")
    if (is.null(opts)) opts = list()
    newopts = list(...)
    for (nm in names(newopts))
        opts[[nm]] = newopts[[nm]]
    options(lsmeans = opts)
    invisible(opts)
}

# Utility that returns TRUE if getOption("lsmeans")[[opt]] is TRUE
.lsm.is.true = function(opt) {
    x = getOption("lsmeans")[[opt]]
    if (is.null(x))  FALSE
    else if (is.logical(x))  x
    else FALSE
}


### Utility to change the internal structure of a ref.grid
### Returned ref.grid object has linfct = I and bhat = estimates
### Primary reason to do this is with transform = TRUE, then can 
### work with linear functions of the transformed predictions
regrid = function(object, transform = TRUE) {
    est = .est.se.df(object, do.se = FALSE)
    estble = !(is.na(est[[1]]))
    object@V = vcov(object)[estble, estble, drop=FALSE]
    object@bhat = est[[1]]
    object@linfct = diag(1, length(estble))
    if(all(estble))
        object@nbasis = matrix(NA)
    else
        object@nbasis = object@linfct[, !estble, drop = FALSE]
    if(transform && !is.null(object@misc$tran)) {
        link = attr(est, "link")
        D = diag(link$mu.eta(object@bhat[estble]))
        object@bhat = link$linkinv(object@bhat)
        object@V = D %*% tcrossprod(object@V, D)
        inm = object@misc$inv.lbl
        if (!is.null(inm))
            object@misc$estName = inm
        object@misc$tran = object@misc$inv.lbl = NULL
    }
    # Nix out things that are no longer needed or valid
    object@grid$.offset. = object@misc$offset.mult =
        object@misc$estHook = object@misc$vcovHook = NULL
    object
}


### S4 methods
## use S3 for this setMethod("summary", "ref.grid", summary.ref.grid)
setMethod("show", "ref.grid", function(object) str.ref.grid(object))
