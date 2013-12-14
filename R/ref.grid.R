# Reference grid code

# TO DO:
# For multivariate respons, expand X via kronecker(Iden, X) 
# where Iden is a kxk identity matrix. Also expand the nbasis
# by replication, and create leves for the response "factor" to
# add to the list

# S4 class definition:
setClass("ref.grid", representation (
    predictors = "character",
    responses = "character",
    grid = "data.frame", 
    levels = "list",
    matlevs = "list",
    X = "matrix",
    bhat = "numeric",
    nbasis = "matrix",
    V = "matrix",
    ddfm = "function",
    misc = "list"
))

setMethod("show", "ref.grid", function(object) {
    showlevs = function(x) # internal convenience function
        cat(paste(format(x, digits = 5, justify = "none"), collapse=", "))
    cat("responses: ")
    showlevs(object@responses)
    levs = object@levels
    cat("\npredictors:\n")
    for (nm in object@predictors) {
        cat(paste("    ", nm, " = ", sep = ""))
        if (nm %in% names(object@matlevs)) {
            cat(paste("matrix with constant columns: ", sep=""))
            showlevs(object@matlevs[[nm]])
        }
        else
            showlevs(levs[[nm]])
        cat("\n")
    }
})

# Change to cov.reduce specification: can be...
#     a function: is applied to all covariates
#     named list of functions: applied to those covariates (else mean is used)
#     TRUE - same as mean
#     FALSE - same as function(x) sort(unique(x))

ref.grid <- function(object, at, cov.reduce = mean) {
    # recover the data
    data = .recover.data (object)
    
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
    
    new ("ref.grid", 
         predictors = attr(data, "predictors"), responses = attr(data, "responses"),
         grid = grid, levels = ref.levels, matlevs = matlevs,
         X = basis$X, bhat = basis$bhat, nbasis = basis$nbasis, V = basis$V,
         ddfm = basis$ddfm, misc = basis$misc)
}
