# Reference grid code

# TO DO: Is there a less clunk way to do the mult.levs argument
# in ref.grid? 

# S4 class definition:
setClass("ref.grid", representation (
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
    
    new ("ref.grid",
         roles = list(predictors = attr(data, "predictors"), 
                      responses = attr(data, "responses"), multresp = multresp),
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
    names(result) = c("estimate", "SE", "df")
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

setMethod("summary", "ref.grid", function(object, ...) {
    # for nor...
    .est.se.df(object@linfct, object@bhat, object@nbasis, object@V, object@ddfm, object@misc)
})
