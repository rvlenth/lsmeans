# Reference grid code

# S4 class definition:
setClass("ref.grid", representation(
    grid = "data.frame", 
    levels = "list"
))

setMethod("show", "ref.grid", function(object) {
    levs = object@levels
    cat("\"ref.grid\" object with levels:\n")
    for (nm in names(levs)) {
        cat(paste("\t", nm, " = ", sep = ""))
        cat(paste(levs[[nm]], collapse = ", "))
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
    ref.levels = matdat = list()
    
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
            matdat[[nm]] = if (!is.logical(cov.reduce)) apply(x, 2, cov.reduce, nm)
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
    for (nm in names(matdat))
        grid[[nm]] = matrix(rep(matdat[[nm]], each=nrow(grid)), nrow=nrow(grid))
    
    new ("ref.grid", grid = grid, levels = ref.levels)
}
