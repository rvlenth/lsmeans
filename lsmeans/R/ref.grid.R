# Reference grid code

# TO DO: Is there a less clunky way to do the mult.levs argument
# in ref.grid? 

# S4 class definition:
setClass("ref.grid", representation (
    model.info = "list",
    roles = "list",
    grid = "data.frame", 
    levels = "list",
    matlevs = "list",
    linfct = "matrix",
    bhat = "numeric",
    nbasis = "matrix",
    V = "matrix",
    ddfm = "function",
    misc = "list"   
))
# Note: misc will hold extra params for ddfm, 
# plus at least the following req'd by the summary method
#   estName: column name for the estimate in the summary ["prediction"]
#   *infer: booleans (CIs?, tests?)  [(FALSE,FALSE)]
#   *level: default conf level [.95]
#   *adjust: default adjust method ["none"]
#   famSize: number of means in family
# *starred ones can be provided as arguments to summary

# Change to cov.reduce specification: can be...
#     a function: is applied to all covariates
#     named list of functions: applied to those covariates (else mean is used)
#     TRUE - same as mean
#     FALSE - same as function(x) sort(unique(x))

ref.grid <- function(object, at, cov.reduce = mean, mult.levs) {
    # recover the data
    data = recover.data (object)
    
    # find out if any variables are coerced to factors
    anm = all.names(attr(data, "terms"))    
    coerced = anm[1 + grep("factor|ordered", anm)]
    
    # convenience functions
    sort.unique = function(x) sort(unique(x))
    if(is.logical(cov.reduce)) 
        if(cov.reduce[1]) cov.reduce = mean
        else              cov.reduce = sort.unique
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
        # mentioned in 'at' list
        if (!missing(at) && !is.null(at[[nm]]))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x))
                xlev[[nm]] = ref.levels[[nm]] = levels(x)
        # matrices
        else if (is.matrix(x)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matlevs[[nm]] = apply(x, 2, cr, nm)
            # if cov.reduce returns a vector, average its columns
            if (is.matrix(matlevs[[nm]]))
                matlevs[[nm]] = apply(matlevs[[nm]], 2, mean)
        }
        # covariate not mentioned in 'at'
        else {
            # single numeric pred but coerced to a factor - use unique values
            if (length(grep(nm, coerced)) > 0)             
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
    
    basis = lsm.basis(object, attr(data, "terms"), xlev, grid)
    
    multresp = list()
    ylevs = basis$misc$ylevs
    if(!is.null(ylevs)) { # have a multivariate situation
        if (missing(mult.levs)) {
            yname = multresp = "rep.meas"
            ref.levels[[yname]] = ylevs
        }
        else {
            k = prod(sapply(mult.levs, length))
            if (k != length(ylevs)) 
                stop("supplied 'mult.levs' is of different length than that of multivariate response")
            for (nm in names(mult.levs))
                ref.levels[[nm]] = mult.levs[[nm]]
            multresp = names(mult.levs)
        }
    }
    basis$misc$estName = "prediction"
    basis$misc$infer = c(FALSE,FALSE)
    basis$misc$level = .95
    basis$misc$adjust = "none"
    basis$misc$famSize = nrow(grid)
    
    new ("ref.grid",
         model.info = list(call = attr(data,"call"), terms = attr(data, "terms"), xlev = xlev),
         roles = list(predictors = attr(data, "predictors"), 
                      responses = attr(data, "responses"), 
                      multresp = multresp, xlev = xlev),
         grid = grid, levels = ref.levels, matlevs = matlevs,
         linfct = basis$X, bhat = basis$bhat, nbasis = basis$nbasis, V = basis$V,
         ddfm = basis$ddfm, misc = basis$misc)
}

# utility fcn to get est's, std errors, and df
# returns a data.frame
.est.se.df = function(linfct, bhat, nbasis, V, ddfm, misc) {
    active = which(!is.na(bhat))
    bhat = bhat[active]
    result = apply(linfct, 1, function(x) {
        estble = if(is.na(nbasis[1]))
            TRUE
        else {
            chk = t(nbasis) %*% x
            x = x[active]
            all(abs(chk) <= 1e-6)      
        }
        if (estble) {
            est = sum(bhat * x)
            se = sqrt(sum(x * V %*% x))
            df = ddfm(x, misc)
            c(est, se, df)
        }
        else c(NA,NA,NA)
    })
    result = as.data.frame(t(result))
    names(result) = c(misc$estName, "SE", "df")
    result
}

### =========== Methods for ref.grid class =============================

setMethod("show", "ref.grid", function(object) {
    showlevs = function(x) # internal convenience function
        cat(paste(format(x, digits = 5, justify = "none"), collapse=", "))
    #cat("responses: ")
    #showlevs(object@responses)
    levs = object@levels
    cat("'ref.grid' object with these variables:\n")
    for (nm in union(object@roles$predictors, object@roles$multresp)) {
        cat(paste("    ", nm, " = ", sep = ""))
        if (nm %in% setdiff(names(object@matlevs), object@roles$multresp)) {
            cat("matrix with constant columns: ")
            showlevs(object@matlevs[[nm]])
        }
        else if (nm %in% object@roles$multresp) {
            cat("multivariate response with levels: ")
            showlevs(levs[[nm]])
        }
        else
            showlevs(levs[[nm]])
        cat("\n")
    }
})

# Utility to parse the dots argument
.getPref = function(arg, dots, default) {
    if (is.null(dots[[arg]])) default
    else dots[[arg]]
}

# utility to compute an adjusted p value
.adj.p.value = function(t, df, adjust, fam.size) {
    n.contr = sum(!is.na(t))
    abst = abs(t)
    unadj.p = 2*pt(abst, df, lower.tail=FALSE)
    if (adjust %in% p.adjust.methods)
        p.adjust(unadj.p, adjust, n = n.contr)
    else switch(adjust,
        auto = unadj.p,
        tukey = ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),
        sidak = 1 - (1 - 2*pt(abst, df, lower.tail=FALSE))^n.contr,
        scheffe = pf(t^2/(fam.size-1), fam.size-1, df, lower.tail=FALSE),
        stop("adjust method '", adjust, "' not implemented")
    ) 
}


setMethod("summary", "ref.grid", function(object, ...) {
    # for now...
    result = .est.se.df(object@linfct, object@bhat, object@nbasis, object@V, object@ddfm, object@misc)
    # figure out factors w/ more than one level
    nlev = sapply(object@levels, length)
    lbls = object@grid[which(nlev > 1)]
    zFlag = (all(is.na(result$df)))
    
    # Add any extras
    dots = list(...)
    infer = .getPref("infer", dots, object@misc$infer)
    if(infer[1]) { # add CIs
        level = .getPref("level", dots, object@misc$level)
        quant = 1 - (1 - level)/2
        cv = if(zFlag) qnorm(quant) else qt(quant, result$df)
        cnm = if (zFlag) c("asymp.LCL", "asymp.UCL") else c("lower.CL","upper.CL")
        result[[cnm[1]]] = result[[1]] - cv*result$SE
        result[[cnm[2]]] = result[[1]] + cv*result$SE
    }
    if(infer[2]) { # add tests
        cnm = ifelse (zFlag, "z.ratio", "t.ratio")
        t.ratio = result[[cnm]] = result[[1]] / result$SE
        adj = .getPref("adjust", dots, object@misc$adjust)
        result$p.value = .adj.p.value(t.ratio, result$df, adj, object@misc$famSize)
    }
    summ = cbind(lbls, result)
    by = object@misc$by
    if (!is.null(by)) 
        attr(summ, "by.vars") = by
    class(summ) = c("summary.rg", "data.frame")
    summ
})
