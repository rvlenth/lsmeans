### Code for an enhancement of 'glht' in 'multcomp' package
### Provides for using 'lsm' in similar way to 'mcp'
### This is implemented via the class "lsmlf" -- linear functions for lsmeans
### (also oddly reminiscent of an old Lucky Strike commercial, LSMFT)

# lsm(specs) will be used as 'linfct' argument in glht
# all we need to do is class it and save the arguments
lsm <- function(...) {
    result <- list(...)
    class(result) <- "lsmlf"
    result
}

# here is S3 method for class "lsmlf" in glht
glht.lsmlf <- function(model, linfct, ...) {
    # linfct has the arguments to pass to lsmeans. We just need to fill it out and call
    linfct$object <- model
    linfct$lf <- TRUE
    
    # In this implementation, we'll NOT accommodate lists
    spec <- linfct$specs
    if (is.list(spec)) linfct$specs <- spec[[1]]
    
    # Get the linear function
    lf <- do.call("lsmeans", linfct)
    
    # Just use the last result - will be linfct for lsmeans if no rhs, else for the contrasts
    glht(model, linfct = lf[[length(lf)]], ...)
}


