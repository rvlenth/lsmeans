# Utility to pick out the args that can be passed to a function
.args.for.fcn = function(fcn, args) {
    oknames = names(as.list(args(fcn)))
    mat = pmatch(names(args), oknames)
    args = args[!is.na(mat)]
    mat = mat[!is.na(mat)]
    names(args) = oknames[mat]
    args
}

# Create a list and give it class class.name
.cls.list <- function(class.name, ...) {
    result <- list(...)
    class(result) <- c(class.name, "list")
    result
}

setMethod("show", "lsmobj", function(object) print(summary(object)) )



### lsmeans S3 generics ...
### I am opting to use S3 methods, cascaded for two arguments
### rather than messing with S4 methods

lsmeans = function(object, specs, ...)
    UseMethod("lsmeans", specs)

# 
lsmeans.default = function(object, specs, ...) {
    rgargs = .args.for.fcn(ref.grid, list(object=object, ...))
    RG = do.call("ref.grid", rgargs)
    lsargs = list(object = RG, specs = specs, ...)
    for (nm in names(rgargs)[-1]) lsargs[[nm]] = NULL
    do.call("lsmeans", lsargs)###lsmeans(RG, specs, ...)
}

# signature = (ALL, ALL) <==> (some model object, ALL)
# Look for args: at, cov.reduce = mean, mult.levs
# lsmeans = function(object, specs, ...) {
#     rgargs = .args.for.fcn(ref.grid, list(object=object, ...))
#     RG = do.call("ref.grid", rgargs)
#     lsargs = list(object = RG, specs = specs, ...)
#     for (nm in names(rgargs)[-1]) lsargs[[nm]] = NULL
#     do.call("lsmeans", lsargs)###lsmeans(RG, specs, ...)
# }
# setGeneric("lsmeans")

#setMethod("lsmeans", signature(object="ANY", specs="formula"),
lsmeans.formula =
function(object, specs, contr.list, trend, ...) {
    if (!missing(trend))
        return(lstrends(object, specs, var=trend, ...))
    
    if(length(specs) == 2) { # just a rhs
        by = .find.by(as.character(specs[2]))
        lsmeans(object, all.vars(specs), by = by, ...)
    }
    else {
#        lsms = lsmeans(object, all.vars(specs[-2]), ...)
        contr.spec = all.vars(specs[-3])[1]
        by = .find.by(as.character(specs[3]))
        # Handle old-style case where contr is a list of lists
        if (!missing(contr.list)) {
            cmat = contr.list[[contr.spec]]
            if (!is.null(cmat))
                contr.spec = cmat
        }
        lsmeans(object, specs = all.vars(specs[-2]), 
                by = by, contr = contr.spec, ...)
    }
}

# List of specs
#setMethod("lsmeans", signature(object="ANY", specs="list"),
lsmeans.list = function(object, specs, ...) {
    result = list()
    nms = names(specs)
    # Format a string describing the results
    .make.desc = function(meth, pri, by) {
        pri = paste(pri, collapse = ", ")
        desc = paste(meth, "of", pri)
        if (!is.null(by)) {
            by = paste(by, collapse = ", ")
            desc = paste(desc, "|", by)
        }
        desc
    }
    
    for (i in seq_len(length(specs))) {
        res = lsmeans(object=object, specs = specs[[i]], ...)
        nm = nms[i]
        if (is.data.frame(res)) { # happens e.g. when cld is used
            if (is.null(nm))
                nm = .make.desc("summary", attr(res, "pri.vars"), attr(res, "by.vars"))
            result[[nm]] = res
        }
        else if (is.list(res)) {
            for (j in seq_len(length(res))) {
                m = res[[j]]@misc
                if (is.null(nm))
                    names(res)[j] = .make.desc(m$methDesc, m$pri.vars, m$by.vars)
                else
                    names(res)[j] = paste(nm, m$methDesc)
            }
            result = c(result,res)
        }
        else{
            if (is.null(nm))
                nm = .make.desc(res@misc$methDesc, res@misc$pri.vars, res@misc$by.vars)
            result[[nm]] = res
        }
    }
    class(result) = c("lsm.list", "list")
    result
}


# Generic for after we've gotten specs in character form
lsmeans.character = function(object, specs, ...)
    UseMethod("lsmeans.character", object)

# Needed for model objects
lsmeans.character.default = function(object, specs, ...)
    lsmeans.default(object, specs, ...)

# Method for a ref.grid -- all methods will get us here eventually
#setMethod("lsmeans", signature(object="ref.grid", specs="character"), 
lsmeans.character.ref.grid = function(object, specs, by = NULL, 
         fac.reduce = function(coefs) apply(coefs, 2, mean), contr, ...) {
    
    RG = object
    facs = union(specs, by)
    
    # Figure out the structure of the grid
    dims = sapply(RG@levels, length)
    row.idx = array(seq_len(nrow(RG@linfct)), dims)
    
    # Get the required factor combs
    levs = list()
    for (f in facs) {
        levs[[f]] = RG@levels[[f]]
        if (is.null(levs[[f]]))
            stop(paste("No variable named", f, "in the reference grid"))
    }
    combs = do.call("expand.grid", levs)
    K = plyr::alply(row.idx, match(facs, names(RG@levels)), function(idx) {
        fac.reduce(RG@linfct[idx, , drop=FALSE])
    })
        
    linfct = t(as.matrix(as.data.frame(K)))
    row.names(linfct) = NULL
    
    if(.some.term.contains(union(facs, RG@roles$trend), RG@model.info$terms))
        message("NOTE: Results may be misleading due to involvement in interactions")
    
    # Figure offset, if any
    if (".offset." %in% names(RG@grid)) {
        offset = plyr::alply(row.idx, match(facs, names(RG@levels)), function(idx) {
            fac.reduce(RG@grid[idx, ".offset.", drop=FALSE])
        })
        combs[[".offset."]] = unlist(offset)
    }
    
    # Figure out which factors have been averaged over
    nlev = sapply(RG@levels, length)
    avgd.over = setdiff(names(nlev[nlev > 1]), facs)
    
    
    RG@roles$responses = character()
    RG@misc$famSize = nrow(linfct)
    RG@misc$estName = "lsmean"
    RG@misc$adjust = "none"
    RG@misc$infer = c(TRUE,FALSE)
    RG@misc$pri.vars = setdiff(facs, by)
    RG@misc$by.vars = by
    RG@misc$avgd.over = union(RG@misc$avgd.over, avgd.over)
    RG@misc$methDesc = "lsmeans"
    RG@roles$predictors = names(levs)
    result = new("lsmobj", RG, linfct = linfct, levels = levs, grid = combs)
    
    if (missing(contr))
        result
    
    else { # return a list with lsmeans and contrasts
        if (is.character(contr) && contr == "cld") {
        # TO DO: provide for passing dots to cld                
            return(cld(result, by = by))
        }
        ctrs = contrast(result, method = contr, by, ...)
        .cls.list("lsm.list", lsmeans = result, contrasts = ctrs)
    }
}


# Summary method for an lsm.list
summary.lsm.list <- function(object, ...)
    lapply(object, function(x) {
        if (inherits(x, "summary.ref.grid"))  x
        else summary(x, ...)
    })

print.lsm.list <- function(x, ...) 
    print(summary(x, ...))


# utility to parse 'by' part of a formula
.find.by = function(rhs) {
    b = strsplit(rhs, "\\|")[[1]]
    if (length(b) > 1) 
        all.vars(as.formula(paste("~",b[2])))
    else NULL
}

### 'contrast' S3 generic and method
contrast = function(object, ...)
    UseMethod("contrast")
              
contrast.ref.grid = function(object, method = "eff", by, adjust, ...) {
    args = object@grid
    args[[".offset."]] = NULL # ignore the offset in labels, etc.
    if(missing(by)) 
        by = object@misc$by.vars
    if (!is.null(by)) {
        by.rows = .find.by.rows(args, by)
        bylevs = args[, by, drop=FALSE]
        args = args[by.rows[[1]], , drop=FALSE]
        for (nm in by) args[[nm]] = NULL
    }
    args$sep = ","
    levs = do.call("paste", args)
    
    
    if (is.list(method)) {
        cmat = as.data.frame(method)
        method = function(levs) cmat
    }
    else if (is.character(method)) {
        fn = paste(method, "lsmc", sep=".")
        method = if (exists(fn, mode="function")) 
            get(fn) 
        else 
            stop(paste("Contrast function '", fn, "' not found", sep=""))
    }
    # case like in old lsmeans, contr = list
    else if (!is.function(method))
        stop("'method' must be a function or the basename of an '.lsmc' function")
    
    # Get the contrasts; this should be a data.frame
    cmat = method(levs, ...)
    if (!is.data.frame(cmat))
        stop("Contrast function must provide a data.frame")
    else if (nrow(cmat) != nrow(args))
        stop("Nonconforming number of contrast coefficients")
    
    if (is.null(by)) {
        linfct = t(cmat) %*% object@linfct
        grid = data.frame(.contrast.=names(cmat))
        if (".offset." %in% names(object@grid))
            grid[[".offset."]] = t(cmat) %*% object@grid[[".offset."]]
    }
    
    # NOTE: The kronecker thing here is nice and efficient but depends
    # on the grid being regular -- same number of rows for each 'by' case
    # If you ever want to expand to irregular grids, this block will
    # have to change, but everything else is probably OK.
    else {
        tcmat = kronecker(diag(rep(1,length(by.rows))), t(cmat))
        linfct = tcmat %*% object@linfct[unlist(by.rows), ]
        tmp = expand.grid(con= names(cmat), by = seq_len(length(by.rows)))###unique(by.id))
        grid = data.frame(.contrast. = tmp$con)
        n.each = ncol(cmat)
        row.1st = sapply(by.rows, function(x) x[1])
        xlevs = list()
        for (v in by)
            xlevs[[v]] = rep(bylevs[row.1st, v], each=n.each)
        grid = cbind(grid, as.data.frame(xlevs))
        if (".offset." %in% names(object@grid))
            grid[[".offset."]] = tcmat %*% object@grid[unlist(by.rows), ".offset."]
    }
    
    # Rename the .contrast. column -- ordinarily to "contrast",
    # but otherwise a unique variation thereof
    n.prev.con = length(grep("^contrast[0-9]?", names(grid)))
    con.col = grep("\\.contrast\\.", names(grid))
    con.name = paste("contrast", 
                     ifelse(n.prev.con == 0, "", n.prev.con), sep="")
    names(grid)[con.col] = con.name
    
    row.names(linfct) = NULL
    misc = object@misc
    misc$estName = "estimate"
    misc$methDesc = attr(cmat, "desc")
    misc$famSize = size=nrow(args)
    if (missing(adjust)) adjust = attr(cmat, "adjust")
    if (is.null(adjust)) adjust = "none"
    misc$adjust = adjust
    misc$infer = c(FALSE, TRUE)
    misc$by.vars = by
    # zap the transformation info except in very special cases
    if (!is.null(misc$tran)) {
        misc$orig.tran = misc$tran
        # anything other than (-1,0,1)?
        non.comp = setdiff(zapsmall(unique(as.matrix(cmat))), c(-1,0,1)) 
        if(length(non.comp) == 0 && (misc$tran %in% c("log", "logit"))) {
            misc$orig.inv.lbl = misc$inv.lbl
            misc$inv.lbl = ifelse(misc$tran == "logit", "odds.ratio", 
                                  paste(misc$inv.lbl,"ratio",sep="."))
            misc$tran = "log"
        }
        else
            misc$tran = NULL
    }
    
    object@roles$predictors = "contrast"
    levels = list()
    for (nm in setdiff(names(grid), ".offset."))
        levels[[nm]] = unique(grid[[nm]])
        
    new("lsmobj", object, linfct=linfct, levels=levels, grid=grid, misc=misc)
}

# return list of row indexes in tbl for each combination of by
# tbl should be a data.frame
.find.by.rows = function(tbl, by) {
    if (any(is.na(match(by, names(tbl)))))
        stop("'by' variables are not all in the grid")    
    bylevs = tbl[ , by, drop=FALSE]
    by.id = do.call("paste", bylevs)
    uids = unique(by.id)
    lapply(uids, function(id) which(by.id == id))
}


# confint method
confint.ref.grid = function(object, parm, level=.95, ...) {
    summary(object, infer=c(TRUE,FALSE), level=level, ...)
}

# test S3 generic and method
test = function(object, parm, ...) {
    UseMethod("test")
}
test.ref.grid = function(object, parm, ...) {
    summary(object, infer=c(FALSE,TRUE), ...)
}

# pairs method
pairs.ref.grid = function(x, ...) {
    object = x # for my sanity
    contrast(object, method = "pairwise", ...)
}


### lstrends function
lstrends = function(model, specs, var, delta.var=.01*rng, data, ...) {
    estName = paste(var, "trend", sep=".") # Do now as I may replace var later

    if (missing(data)) {
        data = try(recover.data (model, data = NULL))
        if (inherits(data, "try-error"))
            stop("Possible remedy: Supply the data used in the 'data' argument")
    }
    else # attach needed attributes to given data
        data = recover.data(model, data = data)
    
    x = data[[var]]
    fcn = NULL   # differential
    if (is.null(x)) {
        fcn = var
        var = all.vars(as.formula(paste("~",var)))
        if (length(var) > 1)
            stop("Can only support a function of one variable")
        else {
            x = data[[var]]
            if (is.null(x)) stop("Variable '", var, "' is not in the dataset")            
        }
    }
    rng = diff(range(x))
    if (delta.var == 0)  stop("Provide a nonzero value of 'delta.var'")
    
    RG = ref.grid(model, ...)
    grid = RG@grid
    grid[[var]] = grid[[var]] + delta.var
    
    basis = lsm.basis(model, attr(data, "terms"), RG@roles$xlev, grid)
    if (is.null(fcn))
        newlf = (basis$X - RG@linfct) / delta.var
    else {
        y0 = with(RG@grid, eval(parse(text = fcn)))
        yh = with(grid, eval(parse(text = fcn)))
        diffl = (yh - y0)
        if (any(diffl == 0)) warning("Some differentials are zero")
        newlf = (basis$X - RG@linfct) / diffl
    }
    
    RG@linfct = newlf
    RG@roles$trend = var
    args = list(object=RG, specs=specs, ...)
    args$at = args$cov.reduce = args$mult.levs = NULL
    result = do.call("lsmeans", args)
    if (is.list(result)) {
        names(result)[1] = "lstrends"
        if (is(result[[1]], "ref.grid")) {
            result[[1]]@misc$estName = estName
            result[[1]]@misc$methDesc = "trends"
        }
    }
    else {
        result@misc$estName = estName
        result@misc$methDesc = "trends"
    }
    
    # No transformation info here
    if (!is.null(result@misc$tran)) {
        result@misc$orig.tran = result@misc$tran
        result@misc$tran = NULL
    }
    
    result
}


# Check if model contains a term containing all elts of facs
# Note: if an lstrends call, we want to include trend var in facs
# terms is terms() component of model
.some.term.contains = function(facs, terms) {
    for (trm in attr(terms, "term.labels")) {
        if(all(sapply(facs, function(f) length(grep(f,trm))>0)))
            if (length(all.vars(as.formula(paste("~",trm)))) > length(facs)) 
                return(TRUE)
    }
    return(FALSE)
}















#===================================================================
#===================================================================
#===================================================================
#===================================================================
#===================================================================
#===================================================================
#===================================================================
#===================================================================
#===================================================================
### Old version of lsmeans
.old.lsmeans = function(object, specs, adjust=c("auto","tukey","sidak","scheffe",p.adjust.methods), conf = .95, 
                   at, trend, contr=list(),
                   cov.reduce = function(x, name) mean(x), 
                   fac.reduce = function(coefs, lev) apply(coefs, 2, mean), 
                   glhargs = NULL, lf = FALSE, mlf = rep(1, nresp) / nresp,
                   ...) 
{
    
    if (missing(specs)) stop("Must specify specs, e.g. 'pairwise ~ treatment'")
    if(!is.null(glhargs)) { # we'll pass contrasts to glht, if multcomp installed; else don't, and warn
        # Version 1.10-2 start requiring multcomp        
        #         if(!require("multcomp")) {
        #             glhargs = NULL
        #             warning("'glhargs' option disabled because 'multcomp' package not installed")
        #         }
        # force integer df since mvtnorm no longer supports fractional
        #        else {
        # I choose to round up if it's within .2 of next integer
        if (!is.null(glhargs$df)) glhargs$df = as.integer(max(1, .2 + glhargs$df))
        #        }
    }
    
# added 6-18-2013 - allow cov.reduce to be logical
    if (is.logical(cov.reduce)) {
        if (cov.reduce) cov.reduce = function(x, name) mean(x)
        else cov.reduce = function(x, name) sort(unique(x))
    }
    
    # for later use
    adjtbl = c("auto","tukey","sidak","scheffe",p.adjust.methods)
    no.adj = pmatch("none", adjtbl)
    adj = pmatch(adjust, adjtbl)[1]
    if (is.na(adj)) {
        adj = 1
        warning("Unknown or non-unique `adjust' method -- automatic method will be used")
    }
    autoadj = (adj == 1)
    
    trend.flag = !missing(trend)
    
# Get model formula and index of response
    if (inherits(object, "gls")) {
        Terms = getCovariateFormula(object)
        yidx = 0
    }
    else {
        Terms = terms(object)
        yidx = attr(Terms, "response")
    }
### was    Terms = delete.response(terms(object))
### but we need to keep track of complete cases, including y values
### (Correction 2-13-13, version 1.06-05)    
    
    # get the pure formula w/o extra stuff
    formrhs = formula(Terms)
    
# ddfm will be replaced with a function of k if there is a way to get denom df    
    ddfm = adjV = NULL
    
# Figure out thecall (fixed effects part of model), bhat (=coefs), contrasts attr
    if (inherits(object, "mer") || inherits(object, "merMod")) {
        if(!isLMM(object) && !isGLMM(object)) 
            stop("Can't handle a nonlinear mixed model")
        thecall = slot(object, "call")
        bhat = fixef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (isLMM(object)) {
            if (require("pbkrtest")) {
                adjV = vcovAdj(object, 0)
                ddfm = function(k) .KRdf.mer (adjV, V, k)
            }
            else warning("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
        }
    }
    else if (inherits(object, "lme")) {
        thecall = object$call
        bhat = fixef(object)
        contrasts = object$contrasts
    }
    else if (inherits(object, "gls")) {
        thecall = object$call
        bhat = coef(object)
        contrasts = object$contrasts
        the.df = object$dims$N - object$dims$p
        ddfm = function(k) the.df
    }
    else if (inherits(object, "lm")) {  ## Also OK for aov, glm, rlm (MASS). Not lqs (but close?)
        thecall = object$call
        bhat = coef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (!(family(object)$family %in% c("binomial", "poisson")))
            if (!is.na(object$df.residual)) 
                ddfm = function(k) object$df.residual
    }
    else
        stop(paste("Can't handle an object of class", class(object)[1]))
    
    
    # Fixed-effects covariance matrix -- Happily, vcov works the same way for lm, lme, lmer
    if(is.null(adjV)) V = vcov(object)
    else V = adjV
    
    # If a multivariate model, reduce to univariate model y %*% mlf
    if (is.matrix(bhat)) {
        nresp = ncol(bhat)
        Iden = diag(rep(1, nrow(bhat)))
        K = kronecker(matrix(mlf, nrow=1), Iden)
        V = K %*% V %*% t(K)
        bhat = bhat %*% mlf
    }
    
    
    # We'll work only with the non-NA elements of bhat
    used = which(!is.na(bhat))
    not.used = which(is.na(bhat))
    bhat = bhat[used]
    # should match vcov...
    if (length(bhat) != nrow(V)) stop("Something's wrong -- Mismatch between vcov() and non-missing coef() results")
    
    # get basis for non-estimable fcns. If NULL, everything is estimable
    null.basis = NULL
    # currently, below can only happen for lm objects. May need to revisit this
    # if other model objects can produce rank-deficient fits
    if (length(not.used) > 0) {
        # null space of X is same as null space of R in QR decomp
        tR = t(qr.R(object$qr))
        if (ncol(tR) < nrow(tR)) # add columns if not square
            tR = cbind(tR, matrix(0, nrow=nrow(tR), ncol=nrow(tR)-ncol(tR)))
        rank = object$qr$rank
        # last few rows are zero -- add a diagonal
        for (i in (rank+1):nrow(tR)) 
            tR[i,i] = 1
        null.basis = qr.resid(qr(tR[, seq_len(rank)]), tR[, -seq_len(rank)])
        if (!is.matrix(null.basis)) null.basis = matrix(null.basis, ncol=1)
        # permute the rows via pivot
        null.basis[object$qr$pivot, ] = null.basis
    }
    
    # All the variables in the model
    nm = all.vars(formrhs)

# Figure out if any are coerced to factor or ordered
    anm = all.names(formrhs)    
    coerced = anm[1 + grep("factor|ordered", anm)]
    
# Obtain a simplified formula -- needed to recover the data in the model    
    form = as.formula(paste("~", paste(nm, collapse = "+")))
    envir = attr(Terms, ".Environment")
    X = model.frame(form, eval(thecall$data, envir=envir), 
                    subset = eval(thecall$subset, enclos=envir),
                    na.action = na.omit, drop.unused.levels = TRUE)
### Correction 2-13-13, version 1.06-05 -- 
### added drop.unused terms and na.omit to above
    
    # Now X contains the data used to fit the model, w/o any expansions (e.g. poly() calls)

# Start accumulating info for the vars. 
# baselevs has the levels of all factors, or the "at" values for all covariates
# xlev has the factor levels only, for use in model.frame calls
    baselevs = xlev = matdat = list()
    # allow a vector of character strings
    if (is.character(specs)) specs = as.list(specs)
    # allow a single formula
    if (!is.list(specs)) specs = list(specs)
    
    all.var.names = names(X)
 
    for (xname in all.var.names) {
        obj = X[[xname]]
        if (is.factor(obj)) {            
            xlev[[xname]] = levels(obj)
            if (!missing(at) && !is.null(at[[xname]]))
                baselevs[[xname]] = at[[xname]]
            else
                baselevs[[xname]] = levels(obj)
        }
        else if (is.matrix(obj)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matdat[[xname]] = apply(obj, 2, cov.reduce, xname)
        }
        else {
            # single numeric pred but coerced to a factor - use unique values
            if (length(grep(xname, coerced)) > 0)             
                 baselevs[[xname]] = sort(unique(obj))
                
            # Ordinary covariates - summarize if not in 'at' arg
            else {
                if (!missing(at) && !is.null(at[[xname]]))
                    baselevs[[xname]] = at[[xname]]
                else 
                    baselevs[[xname]] = cov.reduce(obj, xname)
            # Keep track of covariates with more than one level
            #    if (length(baselevs[[xname]]) > 1)
            #        mult.covar = c(mult.covar, xname)
                
            }
        }
    }
    
    # If 'trend' present, we'll set up a difference quotient based on 
    # a fraction of the range of the variable
    if (trend.flag) {
        # trend could be a variable or a term; 
        # set up trend.xnm to always be variable names(s)
        if (! is.character(trend))
            stop("'trend' must be of character type")
        trend.xnm = trend[1]
        trend.x = X[[trend.xnm]]
        if (is.null(trend.x)) { 
# 'trend' not found, so hold out for possibility that it is the name of a term
            trend.xnm = try(all.vars(as.formula(paste("~",trend))), silent=TRUE)
            if (inherits(trend.xnm, "try-error")) 
                trend.xnm = character(0)
            trend.h = -1 # flag that trend is a term
        }
        else {
            if (!is.numeric(trend.x))
                stop("'trend' must refer to a numeric predictor")
            trend.h = diff(range(trend.x))*.001
            baselevs[[trend]] = baselevs[[trend]][1] + c(-1,1) * trend.h / 2
            # rearrange ordering so trend variable is first - makes bookkeeping easier
            sidx = match(trend, names(baselevs))
            baselevs = c(baselevs[sidx], baselevs[-sidx])
        }
    }
    
    # Keep the response variable from enlarging the grid, no matter what
    if (yidx > 0) {
        yname = as.character(attr(Terms, "variables")[[1 + yidx]])
        if (!is.na(match(yname, names(baselevs))[1])) 
            baselevs[[yname]] = NA
    }
    
    # OK. Now make a grid of the factor levels of interest, along w/ covariate "at" values
    grid = do.call(expand.grid, baselevs)
    
    # add any matrices
    for (nm in names(matdat))
        grid[[nm]] = matrix(rep(matdat[[nm]], each=nrow(grid)), nrow=nrow(grid))

    # It turns out that numerics coerced to factors are a real pain in the butt when it comes
    # to matching levels. Life will be simpler if we turn them into factors in the X matrix 
    # and update the base levels accordingly with the same labels
    #
    #--- Version 1.10 - I don't think I need this anymore with new matching routine
#     for (var in coerced) {
#         X[[var]] = factor(X[[var]])
#         baselevs[[var]] = levels(X[[var]])
#     }

    # Now make a new dataset with just the factor combs and covariate values we want for prediction
    # WARNING -- This will overwrite X, so get anything you need from X BEFORE we get here
    m = model.frame(Terms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(Terms, m, contrasts.arg = contrasts)

    # Compute the derivatives when 'trend' is provided
    if(trend.flag) {
        if (trend.h < 0) { # 'trend' is a term - fix the coefs
            # wipe-out spaces to make it less fussy
            term.nm = gsub(" ", "", dimnames(X)[[2]])
            trend = gsub(" ", "", trend)
            trend.idx = match(trend, term.nm)
            if (is.na(trend.idx))
                stop(paste("trend value '", trend, "' is neither a variable nor a model term", sep=""))
            prev.X = X
# Symbolic differentiation... Replace columns of X with:
#  1 if it is col for trend
#  0 if it does not contain trend
#  previous version of X[,j] where j is interaction of other predictors            
            for (i in seq_len(ncol(X))) {
                trm = term.nm[i]
                trm.pieces = strsplit(trm, ":")[[1]]
                trm.mat = match(trend, trm.pieces)
                if (is.na(trm.mat))
                    X[, i] = 0
                else {
                    trm.otr = paste(trm.pieces[-trm.mat], collapse=":")
                    trm.ref = match(trm.otr, term.nm)
                    if (is.na(trm.ref))     X[, i] = 0 + (i == trend.idx)
                    else                    X[, i] = prev.X[, trm.ref]
                }
            }
        }
        else { # 'trend' is a variable - do the diff quotient
            evens = 2 * (seq_len(nrow(X)/2))
            X = (X[evens, ] - X[evens-1, ]) / trend.h
            # restore the base levels
            baselevs[[trend]] = baselevs[[trend]][1] + trend.h/2
            grid = grid[evens, , drop=FALSE]
        }
    }
    
    # If necessary revise grid with coerced numeric factors replaced with factor levels
    if (length(coerced) > 0) grid = do.call("expand.grid", baselevs)
    
    # All factors (excluding covariates)
    # version 1.10 - no longer excluding covariates
    allFacs = all.var.names
    
    ### Array of indexes for rows of X, organized by dimensions
    row.indexes = array(seq_len(nrow(X)), sapply(baselevs, length))
    
    
    
    # Get a vector of terms in the model, for checking
    mod.terms = attr(Terms, "term.labels")
        ### was strsplit(as.character(formrhs[2])[[1]], "\\+")[[1]]
    
##### routine returns TRUE iff all elements of facs are contained in a model term with another predictor
    some.term.contains = function(facs) {
        # When checking trend, add that variable to the mix
        if(trend.flag) 
            facs = union(facs, trend.xnm)
        for (trm in mod.terms) {
            flag = all(sapply(facs, function(f) length(grep(f,trm))>0))
            if (flag) 
                if (length(all.vars(as.formula(paste("~",trm)))) > length(facs)) 
                    return(TRUE)
        }
        return(FALSE)
    }
    
    
    # Initialize a list to hold the results to return
    results = list()
    for (i in seq_len(length(specs))) {
        form = specs[[i]]
        # convert a string to a formula
        if (is.character(form)) form = as.formula(paste("~",form))
        if (!inherits(form, "formula"))
            stop(paste("Incorrect formula specification:", form))
        method = byfacs = NULL
        if (length(form) == 3) { # no lhs
            method = all.vars(form[[2]])[1]
            form = form[-2]
        }
        
        # These are the variables involved; and the label to use in the results
        facs = all.vars(form)
        facs.lbl = paste(facs, collapse=":")
        if (some.term.contains(facs)) 
            warning(paste("lsmeans of",facs.lbl,
                          "may be misleading due to interaction with other predictor(s)"))
        
        ln = if (any(sapply(facs, function(nm) length(grep(nm, allFacs)) == 0)))
            stop(paste("Unknown factor(s) in specification:", paste(form, collapse=" ")))
        
        # identify "by" factors (after "|" in formula)
        b = strsplit(as.character(form[2]), "\\|")[[1]]
        if (length(b) > 1) byfacs = all.vars(as.formula(paste("~",b[2])))
        
        
        # create the grid of factor combinations
        levs = list()
        for (f in facs) levs[[f]] = baselevs[[f]]
        combs = do.call("expand.grid", levs)

### New (version 1.10) more efficient derivation of K matrix
#         RI = plyr:::splitter_a(row.indexes, match(facs, names(baselevs)))
#     # Each entry of RI has the row indexes of X
#     # for each combination of facs (in expand.grid order)
#         K = sapply(RI, function(idx) {
#             fac.reduce(X[idx, , drop=FALSE], "")
#         })
# Yet newer version, more efficient, uses an exported plyr function too
        K = plyr::alply(row.indexes, match(facs, names(baselevs)), function(idx) {
                    fac.reduce(X[idx, , drop=FALSE], "")
        })
        K = as.matrix(as.data.frame(K))
        
#--- above code replaces pre-1.10 code below...
#         # For each comb, find the needed lin. comb. of bhat to estimate
#         # (These will end up being the COLUMNS of K)
#         K = apply(combs, 1, function(lev) {
#             matches = apply(grid, 1, function(row) {
#                 if (is.numeric(lev)) 
#                     all(zapsmall(as.numeric(row[facs]) - lev) == 0) 
#                 else
#                     all(row[facs] == lev)
#             })
#             nmat = sum(matches)
#             if (nmat == 0) stop(paste("Can't predict at level", lev, "of", "facs.lbl"))
#             else fac.reduce(X[matches, , drop=FALSE], lev)
#         })
        rnames = dimnames(K)[[2]] = apply(combs, 1, paste, collapse=", ")
        
    #### Here is the fcn I'll call to table an estimate of k'beta
        do.est = function(k) {
            est = se = df = NA
            estimable = TRUE
            if (!is.null(null.basis)) {
                estimable = all(abs(apply(null.basis, 2, function(x) sum(k*x))) < 1e-4)      
            }
            if (estimable) {
                k = k[used]
                est = sum(k * bhat)
                se = sqrt(sum(k * (V %*% k)))
                if (!is.null(ddfm)) df = ddfm(k)
            }
            c(estimate=est, SE=se, df=df)
        }
    
    ##### Compute adjusted p value
    # I added zapsmall in ptukey call, seems to help with its flaky behavior
        adj.p.value = function(t, df, meth, fam.size, n.contr) {
            abst = abs(t)
            if (meth <= 4)
                switch(meth,
                       NA,                                                     # should not happen
                       ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),   # tukey
                       1 - (1 - 2*pt(abst, df, lower.tail=FALSE))^n.contr,     # sidak
                       pf(t^2/(fam.size-1), fam.size-1, df, lower.tail=FALSE)  # scheffe
                )
            else
                p.adjust(2*pt(abst, df, lower.tail=FALSE), adjtbl[meth], n=n.contr)
        }
        
        # LS means
        if (!trend.flag) {
            effname = "lsmean"
            lsmentry = paste(facs.lbl, "lsmeans")
        }
        else {
            effname = paste(trend, "trend", sep=".")
            lsmentry = paste(effname, "by", facs.lbl)
        }

        if (lf) {
            results[[lsmentry]] = t(K)
        }
        else {
            lsms = as.data.frame(t(apply(K, 2, do.est)))
            # fix-up names and get CIs
            names(lsms)[1] = effname
            # include factor levels
            lsms = cbind(combs, lsms)
            if (conf > 1) conf = conf/100 # pct --> frac
            if ((conf < 1) && (conf > .01)) {
                if (is.null(ddfm)) {
                    me = qnorm((1-conf)/2, lower.tail=FALSE) * lsms$SE
                    lsms$asymp.LCL = lsms[[effname]] - me
                    lsms$asymp.UCL = lsms[[effname]] + me
                }
                else {
                    me = qt((1-conf)/2, lsms$df, lower.tail=FALSE) * lsms$SE
                    lsms$lower.CL = lsms[[effname]] - me
                    lsms$upper.CL = lsms[[effname]] + me
                }
            }
            
#             if(sort) {
#                 ord = order(lsms[[effname]])
#                 lsms = lsms[ord, ]
#                 K = K[, ord]
#             }
            
            attr(lsms, "print.row.names") = FALSE
            class(lsms) = c("data.frame.lsm", "data.frame")
            results[[lsmentry]] = lsms
        }        
        
        # Do requested contrasts
        if (! is.null(method)) {
            cld.flag = (method == "cld")
            if (cld.flag)
                method = "pairwise"
            
            # look for contrast fcn
            fn = paste(method, "lsmc", sep=".")
            confcn = if (exists(fn, mode="function")) get(fn) 
                else NULL
            
            # bylist will be a list of subsets of the combs to be contrasted
            if (is.null(byfacs)) bylist = list(seq_len(nrow(combs))) # all in one set
            else {
                bg = list()
                for (f in byfacs) bg[[f]] = baselevs[[f]]
                bygrid = do.call("expand.grid", bg)
                bylist = lapply(seq_len(nrow(bygrid)), function(row) {
                    bylevs = bygrid[row,]
                    if (length(byfacs)>1) flags = apply(combs[ , byfacs], 1, function(r) all(r==bylevs))
                    else flags = combs[,byfacs] == bylevs
                    which(flags)
                })
                bylabs = apply(bygrid, 1, paste, collapse=",")
                bycols = sapply(byfacs, grep, names(combs))
                rnames = combs[ ,-bycols]
                if (!is.null(ncol(rnames))) rnames = apply(rnames, 1, paste, collapse=",")
            }
            
            # list to hold results
            Clist = list()
            zer = rep(0, ncol(K)) #### replaced nrow(lsms)) when lf arg added 
            
            # OK, let's go thru the bylist
            nby = length(bylist)
            for (i in seq_len(nby)) {
                rows = bylist[[i]]
                cl = if(is.null(confcn)) contr[[method]] 
                    else confcn(rnames[rows] , ...)
                if (is.null(cl)) stop(paste("Unknown contrast family:", method))
                clx = lapply(cl, function(cc) {
                    if(length(cc) != length(rows))
                        stop(paste(length(cc), " contrast coefficients in '", method, 
                                   "' when ", length(rows), " were expected", sep=""))
                    ccc = zer; 
                    ccc[rows]=cc; 
                    ccc
                })
                if (nby > 1) names(clx) = paste(names(clx), "|", bylabs[i])
                
                # Currently, we're combining all these sets of contrasts in one table
                # Maybe we want to reconsider? If so, ship most of code below inside the
                # loop, and modify labels
                Clist = c(Clist, clx)
            }
            
            # Make a good label for the contrast table
            methdesc = attr(cl, "desc")
            
            # Try to at least explain why glht screws up
            if (!is.null(null.basis) && !is.null(glhargs)) {
                #glhargs = NULL
                warning("Error may occur in 'glht' due to rank deficiency")
            }
            if (lf || !is.null(glhargs)) { # create linear fcn for glht
                KK = t(sapply(Clist, function(con) {
                    nz = which(abs(con) > .0001)
                    K[ , nz] %*% con[nz, drop = FALSE]    
                }))
                if (lf) {
                    dimnames(KK)[[2]] = row.names(K)
                    ctbl = KK
                }
                else {
                    # If glht gets fixed for rank deficiency, may want to consider checking rows of KK
                    # for estimability (see code in do.est())
                    args = c(list(model=object, linfct=KK[ , used]), glhargs)
                    ctbl = summary(do.call("glht", args))
                }
            }
            else { # internal way of doing contrasts
                if (is.null(methdesc)) methdesc = method

                # Figure out the multiplicity adjustment
                adjattr = attr(cl, "adjust")
                if (autoadj) adj = ifelse(is.null(adjattr), no.adj, pmatch(adjattr, adjtbl))
                if (is.na(adj)) adj = no.adj
                ctbl = as.data.frame(t(sapply(Clist, function(con) {                    
                    nz = which(abs(con) > .0001)
                    k = K[ , nz, drop = FALSE] %*% con[nz]
                    do.est(k)
                })))
                
                # factors for mult adjustments...
                n.fam = nrow(lsms) / nby  ########### WAS sum(!is.na(lsms$lsmean)) / nby               
                n.contr = sum(!is.na(ctbl$estimate))
                if (!is.null(ddfm)) {
                    ctbl$t.ratio = round(ctbl$estimate / ctbl$SE, 5)
                    ctbl$p.value = round(adj.p.value(ctbl$t.ratio, ctbl$df, adj, n.fam, n.contr), 5)
                }
                else {
                    ctbl$z.ratio = round(ctbl$estimate / ctbl$SE, 5)
                    ctbl$p.value = round(adj.p.value(ctbl$z.ratio, 10000, adj, n.fam, n.contr), 5)
                }
                attr(ctbl, "mesg") = if(adj == 2)
                    paste("p values are adjusted using the", adjtbl[adj], "method for", n.fam, "means")
                else if (adj < length(adjtbl))
                    paste("p values are adjusted using the", adjtbl[adj], "method for", n.contr, "tests")
                else "p values are not adjusted"
                attr(ctbl, "print.row.names") = TRUE
                class(ctbl) = c("data.frame.lsm", "data.frame")
                
                if(cld.flag) { ### add cld to lsm table
stop(".old.lsmeans does not include 'cld' support")
#                     lsms = cbind(.group = .get.cld(ctbl, 1-conf, nby), lsms)
#                     attr(lsms, "print.row.names") = FALSE
#                     class(lsms) = c("data.frame.lsm", "data.frame")
#                     results[[lsmentry]] = lsms
                }
            }
            results[[paste(facs.lbl,methdesc)]] = ctbl
        }
    }
    if (!lf)
        class(results) = c("lsm","list")
    results
}

#---removed - not really needed, & apparently not more efficient than splitter_a
# # My version of plyr:::splitter_a, but returns a matrix
# # On return, each column has the elements of .array for each combination of .margins
# # Order of columns is same as obtained using expand.grid with the same variables
# .mysplit = function(.array, .margins) {
#     dims = dim(.array)
#     len = length(dims)
#     if (any(.margins > len))
#         stop ("'.margins' exceeds dimensions of '.array'")
#     prm = c(setdiff(seq_len(len), .margins), .margins)
#     matrix(aperm.default(.array, prm, FALSE), ncol = prod(dims[.margins]))
# }

### S3 print method for "lsm" class - only reason we need this now is to support the 'omit' arg
print.lsm = function(x, omit=NULL, ...) {
    for (i in seq_len(length(x))) {
        if (i %in% omit) next
        cat(paste("$`", names(x)[i], "`\n", sep="")) # mimic print method for lists
        print(x[[i]])
        cat("\n")
    }
    invisible(x)
}

### S3 print method for annotated data.frames
print.data.frame.lsm = function(x, ...) {
    print.data.frame(x, row.names=attr(x, "print.row.names"))
    msg = attr(x, "mesg")
    if (!is.null(msg)) 
        for (j in seq_len(length(msg))) cat(paste("   ", msg[j], "\n"))
    invisible(x)
}

