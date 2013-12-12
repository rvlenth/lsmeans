# Recover data routines


# Workhorse for methods
#! Note: Seems like we don't need xlev, and if we use it, we
#! can get warning messages for models like y ~ factor(dose)
#! Moreover, looking at model.frame code, we won't get to 
#! drop.unused.levels part if xlev is NOT null!
.rd.wh <- function(fcall, xlev, trms) {
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- as.name("model.frame")
    fcall$xlev <- NULL ##OMIT## xlev
    fcall$na.action <- na.omit
    vars <- all.vars(trms) # (length will always be >= 2)
    # Put one var on left - keeps out lhs transformations
    form <- as.formula(paste(vars[1], "~", paste(vars[-1], collapse="+")))
    fcall$formula <- update(trms, form)
    env <- environment(trms)
    if (is.null(env)) 
        env <- parent.frame()
    tbl <- eval(fcall, env, parent.frame())
    attr(tbl, "terms") = trms
    attr(tbl, "predictors") = all.vars(delete.response(trms))
    attr(tbl, "responses") = setdiff(vars, attr(tbl, "predictors"))
    tbl
}

### S3 Methods...
# All .recover.data methods will return a data.frame with at least these 
# additional attrs:
#   attr(, "terms")      - terms component of object
#   attr(, "responses")  - names of response variable
#   attr(, "predictors") - names of predictors

# generic
.recover.data = function(object, ...)
    UseMethod(".recover.data")

# default -- used for classes we don't know about
.recover.data.default = function(object, ...)
    stop(paste("Can't handle a model of class", class(object)[1]))

# stats...
.recover.data.lm = function(object) {
    fcall = object$call
    xlev = object$xlevels
    .rd.wh(fcall, xlev, terms(object))
}
# (also aov is extension of lm when no Error() in model)

# MASS
.recover.data.lqs = function(object)
    .recover.data.lm(object)
# (also rlm is extension of lm)

# nlme ...
.recover.data.lme = function(object)
    .recover.data.lm(object)

.recover.data.gls = function(object) {
    fcall = object$call
    xlev = object$xlevels
    .rd.wh(fcall, xlev, getCovariateFormula(object))
}

# lme4 ...
# If it turns out we need xlev after all, we'll need to fix this
.recover.data.merMod = function(object) {
    if(!isLMM(object) && !isGLMM(object)) 
        stop("Can't handle a nonlinear mixed model")
    fcall = object@call
    xlev = NULL ###object$xlevels
    .rd.wh(fcall, xlev, terms(object))
}

# lme4.0 if I can get it and test it
.recover.data.mer = function(object)
    .recover.data.merMod(object)

