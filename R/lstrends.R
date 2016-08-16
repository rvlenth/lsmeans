### Code for lstrends


### lstrends function
lstrends = function(model, specs, var, delta.var=.01*rng, data, type = "link", ...) {
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
        var = .all.vars(as.formula(paste("~",var)))
        if (length(var) > 1)
            stop("Can only support a function of one variable")
        else {
            x = data[[var]]
            if (is.null(x)) stop("Variable '", var, "' is not in the dataset")            
        }
    }
    rng = diff(range(x))
    if (delta.var == 0)  stop("Provide a nonzero value of 'delta.var'")
    
    RG = ref.grid(model, data = data, ...)
    .save.ref.grid(RG)  # save in .Last.ref.grid, if enabled
    
    grid = RG@grid
    if (!is.null(mr <- RG@roles$multresp)) {
        # use the grid value only for the 1st mult resp (no dupes)
        if (length(mr) > 0)
            grid = grid[grid[[mr]] == RG@levels[[mr]][1], ]
    }
    grid[[var]] = grid[[var]] + delta.var
    
    basis = lsm.basis(model, attr(data, "terms"), RG@roles$xlev, grid, ...)
    if (is.null(fcn))
        newlf = (basis$X - RG@linfct) / delta.var
    else {
        y0 = with(RG@grid, eval(parse(text = fcn)))
        yh = with(grid, eval(parse(text = fcn)))
        diffl = (yh - y0)
        if (any(diffl == 0)) warning("Some differentials are zero")
        newlf = (basis$X - RG@linfct) / diffl
    }
    
    
    # args for lsmeans calls
    args = list(object=RG, specs=specs, ...)
    args$at = args$cov.reduce = args$mult.levs = args$vcov. = NULL
    
    type = .validate.type(type)
    # Save corresp lsmeans object if there is a transformation - needed later if summarizing trends w/ type = "response"
    if ((type == "response") && hasName(RG@misc, "tran"))
        lsmean = do.call(lsmeans, args)
    else
        lsmean = NULL
    
    # Now replace linfct w/ difference quotient
    RG@linfct = newlf
    RG@roles$trend = var
    RG@misc$tran = RG@misc$tran.mult = NULL
    RG@misc$initMesg = ifelse(is.null(lsmean),
                              "Trends are based on the transformed scale",
                              "Trends are based on the response scale, after back-transforming")
    
    # Create a possibly object, using results from associated lsmeans object
    .lsmobj = function(obj, lsm) {
        if (is(obj, "ref.grid") && !is.null(lsm)) { # happens only for tran present, type = "resp"
            prd = .est.se.df(lsm, do.se = FALSE)
            lnk = attr(prd, "link")
            deriv = lnk$mu.eta(prd[[1]])
            obj@linfct = diag(deriv) %*% obj@linfct
        }
        obj
    }
    
    args$object = RG
    result = do.call("lsmeans", args)
    
    if (is.list(result)) {
        names(result)[1] = "lstrends"
        if (is(result[[1]], "ref.grid")) {
            result[[1]]@misc$estName = estName
            result[[1]]@misc$estType = "prediction"
            result[[1]]@misc$methDesc = "trends"
            for (i in seq_along(result))
                result[[i]] = .lsmobj(result[[i]], lsmean)
        }
    }
    else {
        result@misc$estName = estName
        result@misc$estType = "prediction"
        result@misc$methDesc = "trends"
        result = .lsmobj(result, lsmean)
    }
    
    result
}

# # predict method
# predict.lstobj = function(object, type, ...) {
#     if (missing(type))
#         type = .get.predict.type(object@misc)
#     else
#         type = .validate.type(type)
#     object@misc$tran = NULL
#     if (type == "response" && (length(object@deriv) > 0))
#         object@linfct = diag(object@deriv) %*% object@linfct
#     
#     predict.ref.grid(object, type = "link", ...)
# }
# 
# 
# 
# summary.lstobj = function(object, type, ...) {
#     if (missing(type))
#         type = .get.predict.type(object@misc)
#     else
#         type = .validate.type(type)
#     object@misc$tran = NULL
#     if (type == "response") {
#         if(length(object@deriv) > 0) {
#             object@linfct = diag(object@deriv) %*% object@linfct
#             annot = "Trends are derived from the back-transformed response"
#         }
#         else
#             annot = "Trends are obtained without back-transforming"
#     }
#     else
#         annot = NULL
#     
#     result = summary.ref.grid(object, ..., type = "link")
#     if(!is.null(annot))
#         attr(result, "mesg") = c(attr(result, "mesg"), annot)
#     result
# }
# 
# 
# # lsmip method
# lsmip.lstobj = function(object, formula, type, ...) {
#     if (missing(type))
#         type = .get.predict.type(object@misc)
#     else
#         type = .validate.type(type)
#     object@misc$tran = NULL
#     if (type == "response" && (length(object@deriv) > 0))
#         object@linfct = diag(object@deriv) %*% object@linfct
#     
#     lsmip.default(object, formula, type = "link", ...)
# }
# 
