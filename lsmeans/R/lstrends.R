### Code for lstrends


### lstrends function
lstrends = function(model, specs, var, delta.var=.01*rng, data, 
                    transform = c("none", "response"), ...) {
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
    
    RG = orig.rg = ref.grid(model, data = data, ...)
    
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
    
    transform = match.arg(transform)
    
    # Now replace linfct w/ difference quotient
    RG@linfct = newlf
    RG@roles$trend = var
    if(hasName(RG@misc, "tran")) {
        tran = RG@misc$tran
        if (is.list(tran)) tran = tran$name
        if (transform == "response") {
            prd = .est.se.df(orig.rg, do.se = FALSE)
            lnk = attr(prd, "link")
            deriv = lnk$mu.eta(prd[[1]])
            RG@linfct = diag(deriv) %*% RG@linfct
            RG@misc$initMesg = paste("Trends are obtained after back-transforming from the", tran, "scale")
        }
        else
            RG@misc$initMesg = paste("Trends are based on the", tran, "(transformed) scale")
    }
   
    RG@misc$tran = RG@misc$tran.mult = NULL
    RG@misc$estName = estName
    
    .save.ref.grid(RG)  # save in .Last.ref.grid, if enabled
    
    # args for lsmeans calls
    args = list(object=RG, specs=specs, ...)
    args$at = args$cov.reduce = args$mult.levs = args$vcov. = args$data = args$trend = NULL
    do.call("lsmeans", args)
}

