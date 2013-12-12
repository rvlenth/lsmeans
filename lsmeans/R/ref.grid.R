# Reference grid code

# S4 class definition:
setClass("ref.grid", representation(
    predictors = "character",
    responses = "character",
    grid = "data.frame", 
    levels = "list",
    matlevs = "list"
))

setMethod("show", "ref.grid", function(object) {
    showlevs = function(x) # internal convenience function
        cat(paste(format(x, digits = 5), collapse=", "))
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

ref.grid <- function(object, at, cov.reduce = function(x, name) mean(x)) {
    # recover the data
    data = .recover.data (object)
    
    # find out if any variables are coerced to factors
    anm = all.names(attr(data, "terms"))    
    coerced = anm[1 + grep("factor|ordered", anm)]
    
    # initialize empty lists
    ref.levels = matpred = list()
    
    for (nm in attr(data, "responses"))
        ref.levels[[nm]] = mean(data[[nm]])
    
    for (nm in attr(data, "predictors")) {
        x = data[[nm]]
        # mentioned in 'at' list
        if (!missing(at) && !is.null(at[[nm]]))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x))
                ref.levels[[nm]] = levels(x)
        # matrices
        else if (is.matrix(x)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matpred[[nm]] = if (!is.logical(cov.reduce)) apply(x, 2, cov.reduce, nm)
                           else apply(x, 2, mean)
        }
        # covariate not mentioned in 'at'
        else {
            # cov.reduce == FALSE or
            # single numeric pred but coerced to a factor - use unique values
            if ((is.logical(cov.reduce) && !cov.reduce) || length(grep(nm, coerced)) > 0)             
                ref.levels[[nm]] = sort(unique(x))
            
            # Ordinary covariates - summarize
            else 
                ref.levels[[nm]] = cov.reduce(x, nm)
        }
    }
    
    # Now create the reference grid
    grid = do.call(expand.grid, ref.levels)
    # add any matrices
    for (nm in names(matpred))
        grid[[nm]] = matrix(rep(matpred[[nm]], each=nrow(grid)), nrow=nrow(grid))
    
    new ("ref.grid", 
         predictors = attr(data, "predictors"), responses = attr(data, "responses"),
         grid = grid, levels = ref.levels, matlevs = matpred)
}
