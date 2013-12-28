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

# S3 method for an lsmobj
glht.lsmobj <- function(model, linfct, ...) {
    args = list(model=model, linfct=linfct@linfct, ...)
    # add a df value if not supplied
    if (is.null(args$df)) {
        df = summary(linfct)$df
        if(any(!is.na(df))) {
            args$df = max(1, as.integer(mean(df, na.rm=TRUE) + .25))
            message("Note: df set to ", args$df)
        }
        
    }
    do.call("glht", args)
}

# New S3 method for lsmlf objects
glht.lsmlf <- function(model, linfct, ...) {
    # Just grab the arguments passed to lsm and call lsmeans
    linfct$object <- ref.grid(model)
    lsmo <- do.call("lsmeans", linfct)
    if (is.list(lsmo)) 
        lsmo = lsmo[[length(lsmo)]]
    # Then call the method for lsmobj
    glht(model, lsmo, ...)
}

