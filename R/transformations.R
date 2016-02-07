# Code to implement transformations my way

### My modification/expansion of stats:make.link()
### Also, if not found, returns make.link("identity") modified with
##     unknown = TRUE, name = link
## In addition, I make all links truly monotone on (-Inf, Inf) in
##     lieu of valideta
.make.link = function(link) {
    if (link %in% c("logit", "probit", "cauchit", "cloglog", "identity", "log"))
        result = stats::make.link(link)
    else result = switch(link,
         sqrt = { tmp = make.link("sqrt") 
             tmp$linkinv = function(eta) pmax(0, eta)^2
             tmp$mu.eta = function(eta) 2*pmax(0, eta)
             tmp },
         `1/mu^2` = { tmp = make.link("1/mu^2") 
             tmp$linkinv = function(eta) 1/sqrt(pmax(0, eta))
             tmp$mu.eta = function(eta) -1/(2*pmax(0, eta)^1.5)
             tmp },
         inverse = { tmp = make.link("inverse") 
             tmp$linkinv = function(eta) 1/pmax(0, eta)
             tmp$mu.eta = function(eta) -1/pmax(0, eta)^2
             tmp },
         `/` = .make.link("inverse"),
         reciprocal = .make.link("inverse"),
         { # default if not included, flags it as unknown
             tmp = stats::make.link("identity")
             tmp$unknown = TRUE
             tmp$name = link
             tmp
         }
    )
    result
}

# Implementation of additional transformations, typically ones with parameters
# Returns a list like stats::make.link, but often with an additional "param" member
# types:
#       glog: log(mu + param)
make.tran = function(type = c("genlog", "power", "boxcox", "arcsin"), param = 1) {
    type = match.arg(type)
    switch(type,
        genlog = list(
            linkfun = function(mu) log(mu + param),
             linkinv = function(eta) exp(eta) - param,
             mu.eta = function(eta) exp(eta),
             valideta = function(eta) all(eta > -param),
             param = param,
             name = paste0("log(mu + ", round(param,3), ")")
        ),
        power = {
            if (param == 0) return(stats::make.link("log"))
            list(
                linkfun = function(mu) mu^param,
                linkinv = function(eta) pmax(eta, 0)^(1/param),
                mu.eta = function(eta) pmax(eta, 0)^(1/param - 1) / param,
                valideta = function(eta) all(eta > 0),
                param = param,
                name = paste0("mu^(", round(param,3), ")")
            )
        },
        boxcox = {
            min.eta = ifelse(param > 0, -1 / param, -Inf)
            if (param == 0) stats::make.link("log")
            else list(
                linkfun = function(mu) (mu^param - 1) / param,
                linkinv = function(eta) (1 + param * pmax(eta, min.eta))^(1/param),
                mu.eta = function(eta) (1 + param * pmax(eta, min.eta))^(1/param - 1),
                valideta = function(eta) all(eta > min.eta),
                param = param,
                name = paste0("Box-Cox (lambda = ", round(param, 3), ")") 
            )
        },
        arcsin = list(
            linkfun = function(mu) asin(sqrt(mu)),
            linkinv = function(eta) sin(pmax(pmin(eta, pi/2), 0))^2,
            mu.eta = function(eta) sin(2*pmax(pmin(eta, pi/2), 0)),
            valideta = function(eta) all(eta <= pi/2) && all(eta >= 0),
            name = "arcsin(sqrt(mu))"
        )
    )
}
