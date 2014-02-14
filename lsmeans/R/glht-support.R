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


# S3 method for an lsmobj
# Note: model is redundant, really, so can be omitted
glht.lsmobj <- function(model, linfct, by, ...) {
    if (!requireNamespace("multcomp"))
        stop(sQuote("glht")," requires ", dQuote("multcomp"), " to be installed")
    object = linfct # so I don't get confused
    if (missing(model)) 
        model = .cls.list("lsmwrap", object = object)
    args = list(model = model, ...)
    # add a df value if not supplied
    if (is.null(args$df)) {
        df = summary(linfct)$df
        if(any(!is.na(df))) {
            args$df = max(1, as.integer(mean(df, na.rm=TRUE) + .25))
            message("Note: df set to ", args$df)
        }
    }
    if (missing(by)) by = object@misc$by.vars
    
    nms = setdiff(names(object@grid), by)
    lf = object@linfct
    dimnames(lf)[[1]] = as.character(interaction(object@grid[, nms], sep=", "))
    
    if (is.null(by)) {
        args$linfct = lf
        return(do.call("glht", args))
    }
    
    # (else...)
    by.rows = .find.by.rows(object@grid, by)
    result = lapply(by.rows, function(r) {
        args$linfct = lf[r, , drop=FALSE]
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

### as. glht -- convert my object to glht object
as.glht <- function(object, ...)
    UseMethod("as.glht")

as.glht.default <- function(object, ...)
    stop("Cannot convert an object of class ", sQuote(class(object)[1]),
         " to a ", sQuote("glht"), " object")

as.glht.lsmobj <- function(object, ...)
    glht( , object, ...)


# S3 modelparm method for lsmwrap (S3 wrapper for an lsmobj - see glht.lsmobj)
modelparm.lsmwrap <- function(model, coef., vcov., df, ...) {
    object = model$object
    estimable = ! sapply(object@bhat, is.na)
    .cls.list("modelparm", coef = object@bhat, vcov = object@V,
              df = df, estimable = estimable)
}

# S3 methods for glht.list
summary.glht.list = function(object, ...)
    lapply(object, summary, ...)





