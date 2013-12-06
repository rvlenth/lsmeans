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


## For reference, this is model.frame method for 'lm' objects
model.frame.lm <- function (formula, ...) 
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
                        0)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        m <- match(c("formula", "data", "subset", "weights", 
                     "na.action", "offset"), names(fcall), 0L)
        fcall <- fcall[c(1L, m)]
        fcall$drop.unused.levels <- TRUE
        fcall[[1L]] <- as.name("model.frame")
        fcall$xlev <- formula$xlevels
        fcall$formula <- terms(formula)
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env)) 
            env <- parent.frame()
        eval(fcall, env, parent.frame())
    }
    else formula$model
}

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
