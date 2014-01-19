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
        
    # Save the original levels of factors, no matter what
    if (is.factor(x))
        xlev[[nm]] = levels(x)
    
    # Now go thru and find reference levels...
        # mentioned in 'at' list but not coerced
        if (!(nm %in% coerced) && !missing(at) && !is.null(at[[nm]]))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x))
            ref.levels[[nm]] = levels(x)
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

    basis = lsm.basis(object, attr(data, "terms"), xlev, grid)
    
# Here's a complication. If a numeric predictor was coerced to a factor, we had to
# include all its levels in the reference grid, even if altered in 'at'
# Moreover, whatever levels are in 'at' must be a subset of the unique values
# So we now need to subset the rows of the grid and linfct based on 'at'
    problems = if (!missing(at)) intersect(coerced, names(at)) 
               else character(0)
    if (length(problems > 0)) {
        incl.flags = rep(TRUE, nrow(grid))
        for (nm in problems) {
            # get only "legel" levels
            at[[nm]] = round(at[[nm]], 3)
            at[[nm]] = at[[nm]][at[[nm]] %in% round(ref.levels[[nm]],3)]
            # Now which of those are left out?
            excl = setdiff(round(ref.levels[[nm]], 3), at[[nm]])
            for (x in excl)
                incl.flags[round(grid[[nm]] - x, 3) == 0] = FALSE
            ref.levels[[nm]] = at[[nm]]
        }
        if (!any(incl.flags))
            stop("Reference grid is empty due to mismatched levels in 'at'")
        grid = grid[incl.flags, , drop=FALSE]
        basis$X = basis$X[incl.flags, , drop=FALSE]
    }

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
                      multresp = multresp),
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
    cat("'ref.grid' object with variables:\n")
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

# Utility to parse the dots argument - I deprecated this already
# .getPref = function(arg, dots, default) {
#     if (is.null(dots[[arg]])) default
#     else dots[[arg]]
# }

# utility to compute an adjusted p value
.adj.p.value = function(t, df, adjust, fam.size) {
# do a pmatch of the adjust method, case insensitive
    adj.meths = c("tukey", "sidak", "scheffe", p.adjust.methods)
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k))
        stop("Adjust method '", adjust, "' is not recognized")
    adjust = adj.meths[k]
    
    n.contr = sum(!is.na(t))
    abst = abs(t)
    unadj.p = 2*pt(abst, df, lower.tail=FALSE)
    if (adjust %in% p.adjust.methods)
        pval = p.adjust(unadj.p, adjust, n = n.contr)
    else pval = switch(adjust,
        tukey = ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),
        sidak = 1 - (1 - 2*pt(abst, df, lower.tail=FALSE))^n.contr,
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

### S4 "summary" method for ref.grid (and lsmobj)
# originally I did not have infer, level, or adjust args
# and looked for them in ... -- bad idea, I think
setMethod("summary", "ref.grid", function(object, infer, level, adjust) {
    result = .est.se.df(object@linfct, object@bhat, object@nbasis, object@V, object@ddfm, object@misc)

# figure out factors w/ more than one level
    nlev = sapply(object@levels, length)
    lbls = object@grid[which(nlev > 1)]
    if (nrow(object@grid) == 1) # but if only one row, use everything
        lbls = object@grid
    zFlag = (all(is.na(result$df)))
    
### implement my 'variable defaults' scheme    
    if(missing(infer)) infer = object@misc$infer
    if(missing(level)) level = object@misc$level
    if(missing(adjust)) adjust = object@misc$adjust

    mesg = NULL
    if(infer[1]) { # add CIs
        quant = 1 - (1 - level)/2
        cv = if(zFlag) qnorm(quant) else qt(quant, result$df)
        cnm = if (zFlag) c("asymp.LCL", "asymp.UCL") else c("lower.CL","upper.CL")
        result[[cnm[1]]] = result[[1]] - cv*result$SE
        result[[cnm[2]]] = result[[1]] + cv*result$SE
        mesg = paste("Confidence level used:", level)
    }
    if(infer[2]) { # add tests
        cnm = ifelse (zFlag, "z.ratio", "t.ratio")
        t.ratio = result[[cnm]] = result[[1]] / result$SE
        apv = .adj.p.value(t.ratio, result$df, adjust, object@misc$famSize)
        adjust = apv$adjust # matched name in case it was abbreviated
        result$p.value = apv$pval
        mesg = c(mesg, apv$mesg)
    }
    summ = cbind(lbls, result)
    by = object@misc$by
    attr(summ, "by.vars") = by
    attr(summ, "mesg") = mesg
    class(summ) = c("summary.ref.grid", "data.frame")
    summ
})

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
        fp = x$p.value = format(round(x$p.value,4), nsmall=4)
        x$p.value[fp=="0.0000"] = "<.0001"
    }
    just = sapply(x.save, function(col) if(is.numeric(col)) "R" else "L")
    xc = as.matrix(format.data.frame(x, digits=digits, na.encode=FALSE))
    m = apply(rbind(just, names(x), xc), 2, function(x) {
        w = max(sapply(x, nchar))
        if (x[1] == "R") format(x[-(1:2)], width = w, justify="right")
        else format(x[-(1:2)], width = w, justify="left")
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
        lbls = do.call(paste, c(x[,by.vars, drop=FALSE], sep=", "))
        for (lb in unique(lbls)) {
            rows = which(lbls==lb)
            levs = paste(by.vars, "=", xc[rows[1], by.vars])
            cat(paste(paste(levs, collapse=", ")), ":\n", sep="")
            print(m[rows, ], ..., quote=quote, right=right)
            cat("\n")
        }
    }
    
    msg = attr(x, "mesg")
    if (!is.null(msg))
        for (j in seq_len(length(msg))) cat(paste(msg[j], "\n"))
    
    invisible(x.save)
}

print.ref.grid = function(x, ...) {
    args.prt = list(...)
    args.sum = list(object=x)
    for (key in c("infer","level","adjust")) {
        args.sum[[key]] = .getPref(key, args.prt, NULL)
        args.prt[[key]] = NULL
    }
    args.prt$x = do.call("summary", args.sum)
    do.call("print", args.prt)
}
