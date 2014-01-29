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
glht.lsmobj <- function(model, linfct, by, ...) {
    object = linfct # so I don't get confused
    args = list(model=model, ...)
    # add a df value if not supplied
    if (is.null(args$df)) {
        df = summary(linfct)$df
        if(any(!is.na(df))) {
            args$df = max(1, as.integer(mean(df, na.rm=TRUE) + .25))
            message("Note: df set to ", args$df)
        }
    }
    if (missing(by)) by = object@misc$by.vars
    
    if (is.null(by)) {
        args$linfct = object@linfct
        return(do.call("glht", args))
    }
    
    # (else...)
    by.rows = .find.by.rows(object@grid, by)
    result = lapply(by.rows, function(r) {
        args$linfct = object@linfct[r, , drop=FALSE]
        do.call("glht", args)
    })
    bylevs = lapply(by, function(byv) unique(object@grid[[byv]]))
    names(bylevs) = by
    bygrid = do.call("expand.grid", bylevs)
    levlbls = lapply(by, function(byv) paste(byv, "=", bygrid[[byv]]))
    levlbls$sep = ", "
    names(result) = do.call("paste", levlbls)
    class(result) = c("glht.list", "list")
    result
}

# S3 methods for glht.list
summary.glht.list = function(object, ...)
    lapply(object, summary, ...)




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

