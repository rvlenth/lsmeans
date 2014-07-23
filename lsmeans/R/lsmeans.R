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
    rgargs$options = NULL  # don't pass options to ref.grid
    RG = do.call("ref.grid", rgargs)
    lsargs = list(object = RG, specs = specs, ...)
    for (nm in names(rgargs)[-1]) lsargs[[nm]] = NULL
    do.call("lsmeans", lsargs)###lsmeans(RG, specs, ...)
}

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
lsmeans.character.ref.grid = function(object, specs, by = NULL, 
         fac.reduce = function(coefs) apply(coefs, 2, mean), 
         contr, options = getOption("lsmeans")$lsmeans, weights, ...) {
    RG = object
    facs = union(specs, by)
    
    # Figure out the structure of the grid
    freq = RG@grid[[".freq."]]
    dims = sapply(RG@levels, length)
    row.idx = array(seq_len(nrow(RG@linfct)), dims)
    use.mars = match(facs, names(RG@levels)) # which margins to use
    avgd.mars = setdiff(seq_along(dims)[dims>1], use.mars) # margins that we average over
    
    # Reconcile weights, if there are any margins left
    if ((length(avgd.mars) > 0) && !missing(weights)) {
        if (is.character(weights)) {
            if (is.null(freq))
                message("Frequency information not available -- deferring to fac.reduce")
            else {
                wopts = c("equal","proportional","outer","cells","invalid")
                weights = switch(wopts[pmatch(weights, wopts, 5)],
                    equal = rep(1, prod(dims[avgd.mars])),
                    proportional = as.numeric(plyr::aaply(row.idx, avgd.mars,
                                                          function(idx) sum(freq[idx]))),
                    outer = {
                        ftbl = plyr::aaply(row.idx, avgd.mars,
                                           function(idx) sum(freq[idx]), .drop = FALSE)
                        w = N = sum(ftbl)
                        for (d in seq_along(dim(ftbl)))
                            w = outer(w, plyr::aaply(ftbl, d, sum) / N)
                        as.numeric(w)
                    },
                    cells = "fq",
                    invalid = stop("Invalid 'weights' option: '", weights, "'")
                )
            }
        }
        if (is.numeric(weights)) {
            wtmat = diag(weights)
            wtsum = sum(weights)
            if (wtsum <= 1e-8) wtsum = NA
            fac.reduce = function(coefs) {
                if (nrow(coefs) != nrow(wtmat))
                    stop("Nonconforming number of weights -- need ", nrow(coefs))
                apply(wtmat %*% coefs, 2, sum) / wtsum
            }
        }
    }
    
    # Get the required factor combs
    levs = list()
    for (f in facs) {
        levs[[f]] = RG@levels[[f]]
        if (is.null(levs[[f]]))
            stop(paste("No variable named", f, "in the reference grid"))
    }
    combs = do.call("expand.grid", levs)
    if (!missing(weights) && (weights == "fq"))
        K = plyr::alply(row.idx, use.mars, function(idx) {
            fq = RG@grid[[".freq."]][idx]
            apply(diag(fq) %*% RG@linfct[idx, , drop=FALSE], 2, sum) / sum(fq)
        })
    else
        K = plyr::alply(row.idx, use.mars, function(idx) {
            fac.reduce(RG@linfct[idx, , drop=FALSE])
        })
        
    linfct = t(as.matrix(as.data.frame(K)))
    row.names(linfct) = NULL
    
    if(.some.term.contains(union(facs, RG@roles$trend), RG@model.info$terms))
        message("NOTE: Results may be misleading due to involvement in interactions")
    
    # Figure offset, if any
    if (".offset." %in% names(RG@grid)) {
        combs[[".offset."]] = as.numeric(plyr::aaply(row.idx, use.mars, function(idx)
            fac.reduce(as.matrix(RG@grid[idx, ".offset.", drop=FALSE]))))
    }
    
    avgd.over = names(RG@levels[avgd.mars])
    
    # Update .freq column of grid, if it exists
    if (!is.null(freq)) {
        combs[[".freq."]] = as.numeric(plyr::aaply(row.idx, use.mars, 
            function(idx) sum(freq[idx])))
    }
    
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
    
    if(!is.null(options)) {
        options$object = result
        result = do.call("update.ref.grid", options)
    }
    
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
              
contrast.ref.grid = function(object, method = "eff", by, adjust, 
        options = getOption("lsmeans")$contrast, ...) {
    args = object@grid
    args[[".offset."]] = NULL 
    args[[".freq."]] = NULL # ignore auxiliary stuff in labels, etc.
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
        
    result = new("lsmobj", object, linfct=linfct, levels=levels, grid=grid, misc=misc)
    if(!is.null(options)) {
        options$object = result
        result = do.call("update.ref.grid", options)
    }
    result
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

