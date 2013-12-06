# Recover data routines

# initial broader attempt - not all that successful
recover.data = function(call) {
    of.interest = match(c("formula","data","subset","na.action"), names(call), 0)
    of.interest = of.interest[of.interest > 0]
    args = as.list(call[of.interest])
    vars = all.vars(args$formula)
    form = as.formula(paste("~", paste(vars, collapse="+")))
    args$formula = update.formula(args$formula, form)
    do.call("model.frame", args)
}



# The following function returns a data.frame with the variables
# from the given object. 
# It also adds attributes "predictors" and "responses"
# This function is adapted from model.frame.lm

# Currently this works correctly for obj of class "lm",
# examples warp.data, warp.log, warp.with; and warp.att in data attached
rd = function(obj) {
    fcall <- obj$call
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- as.name("model.frame")
    fcall$xlev <- obj$xlevels
    # replace formula by list of vars w/out function calls
    trms <- terms(obj)
    vars <- all.vars(trms) # (length will always be >= 2)
    # Put one var on left - keeps out lhs transformations
    form <- as.formula(paste(vars[1], "~", paste(vars[-1], collapse="+")))
    fcall$formula <- update(trms, form)
    env <- environment(obj$terms)
    if (is.null(env)) 
        env <- parent.frame()
    tbl <- eval(fcall, env, parent.frame())
    attr(tbl, "terms") = NULL
    attr(tbl, "predictors") = all.vars(delete.response(trms))
    attr(tbl, "responses") = setdiff(vars, attr(tbl, "predictors"))
    tbl
}
