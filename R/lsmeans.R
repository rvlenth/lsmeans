lsmeans = function(object, specs, adjust=c("auto","tukey","sidak",p.adjust.methods), conf = .95, 
                   at, contr=list(), 
                   cov.reduce = function(x, name) mean(x), 
                   fac.reduce = function(coefs, lev) apply(coefs, 2, mean), 
                   glhargs=NULL, lf = FALSE,
                   ...) 
{
    
    if (missing(specs)) stop("Must specify specs, e.g. 'pairwise ~ treatment'")
    if(!is.null(glhargs)) { # we'll pass contrasts to glht, if multcomp installed; else don't, and warn
        if(!require("multcomp")) {
            glhargs = NULL
            warning("'glhargs' option disabled because 'multcomp' package not installed")
        }
        # force integer df since mvtnorm no longer supports fractional
        # I choose to round up if it's within .2 of next integer
        else {
            if (!is.null(glhargs$df)) glhargs$df = as.integer(max(1, .2 + glhargs$df))
        }
    }
    
    # for later use
    adjtbl = c("auto","tukey","sidak",p.adjust.methods)
    no.adj = pmatch("none", adjtbl)
    adj = pmatch(adjust, adjtbl)[1]
    if (is.na(adj)) {
        adj = 1
        warning("Unknown or non-unique `adjust' method -- automatic method will be used")
    }
    autoadj = (adj == 1)
    
# Get RHS of model formula
    if (inherits(object, "gls"))
        Terms = getCovariateFormula(object)
    else
        Terms = delete.response(terms(object))
    # get the pure formula w/o extra stuff
    formrhs = formula(Terms)
    
# ddfm will be replaced with a function of k and se if there is a way to get denom df    
    ddfm = adjV = NULL
    
# Figure out thecall (fixed effects part of model), bhat (=coefs), contrasts attr
    if (inherits(object, "mer")) {
        if(dim(object@V)[1] > 0) ### gradient matrix is nontrivial only for nonlinear models
            stop("Can't handle an 'nlmer' model")
        thecall = slot(object, "call")
        bhat = fixef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (length(object@muEta) == 0) { # no glm's allowed...
    ### Would really rather use below instead, but pbkrtest won't cooperate w/ any glm right now
    ###   fam = object@call$family
    ###   if (!is.null(fam) && !(fam %in% c("binomial", "poisson")))
            if (require("pbkrtest")) {
                adjV = vcovAdj(object, 0)
                ddfm = function(k, se) .KRdf.mer (adjV, V, k, se*se)
            }
            else warning("Install package 'pbkrtest' to obtain bias corrections and degrees of freedom")
        }
    }
    else if (inherits(object, "lme")) {
        thecall = object$call
        bhat = fixef(object)
        contrasts = object$contrasts
    }
    else if (inherits(object, "gls")) {
        thecall = object$call
        bhat = coef(object)
        contrasts = object$contrasts
        the.df = object$dims$N - object$dims$p
        ddfm = function(k, se) the.df
    }
    else if (inherits(object, "lm")) {  ## Also OK for aov, glm, rlm (MASS). Not lqs (but close?)
        thecall = object$call
        bhat = coef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (!(family(object)$family %in% c("binomial", "poisson")))
            if (!is.na(object$df.residual)) 
                ddfm = function(k, se) object$df.residual
    }
    else
        stop(paste("Can't handle an object of class", class(object)[1]))
    
    
    # Fixed-effects covariance matrix -- Happily, vcov works the same way for lm, lme, lmer
    if(is.null(adjV)) V = vcov(object)
    else V = adjV
    
    # We'll work only with the non-NA elements of bhat
    used = which(!is.na(bhat))
    not.used = which(is.na(bhat))
    bhat = bhat[used]
    # should match vcov...
    if (length(bhat) != nrow(V)) stop("Something's wrong -- Mismatch between vcov() and non-missing coef() results")
    
    # get basis for non-estimable fcns. If NULL, everything is estimable
    null.basis = NULL
    # currently, below can only happen for lm objects. May need to revisit this
    # if other model objects can produce rank-deficient fits
    if (length(not.used) > 0) {
        # null space of X is same as null space of R in QR decomp
        tR = t(qr.R(object$qr))
        rank = object$qr$rank
        # last few rows are zero -- add a diagonal
        for (i in (rank+1):nrow(tR)) tR[i,i] = 1
        null.basis = qr.resid(qr(tR[, 1:rank]), tR[, -(1:rank)])
        if (!is.matrix(null.basis)) null.basis = matrix(null.basis, ncol=1)
        # permute the rows via pivot
        null.basis[object$qr$pivot, ] = null.basis
    }
    
    # All the variables in the model
    nm = all.vars(formrhs)

# Figure out if any are coerced to factor or ordered
    anm = all.names(formrhs)    
    coerced = anm[1 + grep("factor|ordered", anm)]
   
# Obtain a simplified formula -- needed to recover the data in the model    
    form = as.formula(paste("~", paste(nm, collapse = "+")))
    envir = attr(Terms, ".Environment")
    X = model.frame(form, eval(thecall$data, envir=envir), 
                    subset = eval(thecall$subset, enclos=envir))
    # Now X contains the data used to fit the model, w/o any expansions (e.g. poly() calls)

# Start accumulating info for the vars. 
# baselevs has the levels of all factors, or the "at" values for all covariates
# xlev has the factor levels only, for use in model.frame and check.cells calls
    baselevs = xlev = matdat = list()
    # allow a vector of character strings
    if (is.character(specs)) specs = as.list(specs)
    # allow a single formula
    if (!is.list(specs)) specs = list(specs)
 
    for (xname in names(X)) {
        obj = X[[xname]]
        if (is.factor(obj)) {            
            xlev[[xname]] = levels(obj)
            if (!missing(at) && !is.null(at[[xname]]))
                baselevs[[xname]] = at[[xname]]
            else
                baselevs[[xname]] = levels(obj)
        }
        else if (is.matrix(obj)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matdat[[xname]] = apply(obj, 2, cov.reduce, xname)
        }
        else {
            # single numeric pred but coerced to a factor - use unique values
            if (length(grep(xname, coerced)) > 0)             
                 baselevs[[xname]] = sort(unique(obj))
                
            # Ordinary covariates - summarize if not in 'at' arg
            else {
                if (!missing(at) && !is.null(at[[xname]]))
                    baselevs[[xname]] = at[[xname]]
                else 
                    baselevs[[xname]] = cov.reduce(obj, xname)
            }
        }
    }
    # OK. Now make a grid of the factor levels of interest, along w/ covariate "at" values
    grid = do.call("expand.grid", baselevs)
    # add any matrices
    for (nm in names(matdat))
        grid[[nm]] = matrix(rep(matdat[[nm]], each=nrow(grid)), nrow=nrow(grid))

    # It turns out that numerics coerced to factors are a real pain in the butt when it comes
    # to matching levels. Life will be simpler if we turn them into factors in the X matrix 
    # and update the base levels accordingly with the same labels
    for (var in coerced) {
        X[[var]] = factor(X[[var]])
        baselevs[[var]] = levels(X[[var]])
    }

    # Now make a new dataset with just the factor combs and covariate values we want for prediction
    # WARNING -- This will overwrite X, so get anything you need from X BEFORE we get here
    m = model.frame(Terms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(Terms, m, contrasts.arg = contrasts)
    # use only the columns with non-missing regr coefs
    #--- OMITTED--we'll handle this later X = X[ , used]
    
    # If necessary revise grid with corced numeric factors replaced with factor levels
    if (length(coerced) > 0) grid = do.call("expand.grid", baselevs)
    
    # All factors (excluding covariates)
    allFacs = c(names(xlev), coerced)
    
    
    
    # Get a vector of terms in the model, for checking
    mod.terms = attr(Terms, "term.labels")
        ### was strsplit(as.character(formrhs[2])[[1]], "\\+")[[1]]
    
##### routine returns TRUE iff all elements of facs are contained in a model term with another predictor
    some.term.contains = function(facs) {
        for (trm in mod.terms) {
            flag = all(sapply(facs, function(f) length(grep(f,trm))>0))
            if (flag) 
                if (length(all.vars(as.formula(paste("~",trm)))) > length(facs)) 
                    return(TRUE)
        }
        return(FALSE)
    }
    
    
    # Initialize a list to hold the results to return
    results = list()
    for (i in 1:length(specs)) {
        form = specs[[i]]
        # convert a string to a formula
        if (is.character(form)) form = as.formula(paste("~",form))
        if (!inherits(form, "formula"))
            stop(paste("Incorrect formula specification:", form))
        method = byfacs = NULL
        if (length(form) == 3) { # no lhs
            method = all.vars(form[[2]])[1]
            form = form[-2]
        }
        
        # These are the variables involved; and the label to use in the results
        facs = all.vars(form)
        facs.lbl = paste(facs, collapse=":")
        if (some.term.contains(facs)) 
            warning(paste("lsmeans of",facs.lbl,
                          "may be misleading due to interaction with other predictor(s)"))
        
        ln = if (any(sapply(facs, function(nm) length(grep(nm, allFacs)) == 0)))
            stop(paste("Unknown factor(s) in specification:", paste(form, collapse=" ")))
        
        # identify "by" factors (after "|" in formula)
        b = strsplit(as.character(form[2]), "\\|")[[1]]
        if (length(b) > 1) byfacs = all.vars(as.formula(paste("~",b[2])))
        
        # create the grid of factor combinations
        levs = list()
        for (f in facs) levs[[f]] = baselevs[[f]]
        combs = do.call("expand.grid", levs)
        
        # For each comb, find the needed lin. comb. of bhat to estimate
        # (These will end up being the COLUMNS of K)
        K = apply(combs, 1, function(lev) {
            matches = apply(grid, 1, function(row) {
                #### DEL if (is.numeric(lev)) all(abs(as.numeric(row[facs]) - lev) < .001)  else 
                all(row[facs] == lev)
            })
            nmat = sum(matches)
            if (nmat == 0) stop(paste("Can't predict at level", lev, "of", "facs.lbl"))
            else fac.reduce(X[matches, , drop=FALSE], lev)
        })
        rnames = dimnames(K)[[2]] = apply(combs, 1, paste, collapse=", ")
        
    #### Here is the fcn I'll call to table an estimate of k'beta
        do.est = function(k) {
            est = se = df = NA
            estimable = TRUE
            if (!is.null(null.basis)) {
                estimable = all(abs(apply(null.basis, 2, function(x) sum(k*x))) < 1e-4)      
            }
            if (estimable) {
                k = k[used]
                est = sum(k * bhat)
                se = sqrt(sum(k * (V %*% k)))
                if (!is.null(ddfm)) df = ddfm(k, se)
            }
            c(estimate=est, SE=se, df=df)
        }
    
    ##### Compute adjusted p value
        adj.p.value = function(t, df, meth, fam.size, n.contr) {
            abst = abs(t)
            if (meth <= 3)
                switch(meth,
                       NA,                                                   # should not happen
                       ptukey(sqrt(2)*abst, fam.size, df, lower.tail=FALSE), # tukey
                       1 - (1 - 2*pt(abst, df, lower.tail=FALSE))^n.contr,   # sidak
                )
            else
                p.adjust(2*pt(abst, df, lower.tail=FALSE), adjtbl[meth], n=n.contr)
        }
        
        # LS means
        if (lf) {
            results[[paste(facs.lbl, "lsmeans")]] = t(K)
        }
        else {
            lsms = as.data.frame(t(apply(K,2,do.est)))
            # fix-up names and get CIs
            names(lsms)[1] = "lsmean"
            # include factor levels
            lsms = cbind(combs, lsms)
            if (conf > 1) conf = conf/100 # pct --> frac
            if ((conf < 1) && (conf > .01)) {
                if (is.null(ddfm)) {
                    me = qnorm((1-conf)/2, lower.tail=FALSE) * lsms$SE
                    lsms$asymp.LCL = lsms$lsmean - me
                    lsms$asymp.UCL = lsms$lsmean + me
                }
                else {
                    me = qt((1-conf)/2, lsms$df, lower.tail=FALSE) * lsms$SE
                    lsms$lower.CL = lsms$lsmean - me
                    lsms$upper.CL = lsms$lsmean + me
                }
            }
            attr(lsms, "print.row.names") = FALSE
            class(lsms) = c("data.frame.lsm", "data.frame")
            results[[paste(facs.lbl, "lsmeans")]] = lsms
        }        
        
        # Do requested contrasts
        if (! is.null(method)) {
            # look for contrast fcn
            fn = paste(method, "lsmc", sep=".")
            confcn = if (exists(fn, mode="function")) get(fn) 
                else NULL
            
            # bylist will be a list of subsets of the combs to be contrasted
            if (is.null(byfacs)) bylist = list(1:nrow(combs)) # all in one set
            else {
                bg = list()
                for (f in byfacs) bg[[f]] = baselevs[[f]]
                bygrid = do.call("expand.grid", bg)
                bylist = lapply(1:nrow(bygrid), function(row) {
                    bylevs = bygrid[row,]
                    if (length(byfacs)>1) flags = apply(combs[ , byfacs], 1, function(r) all(r==bylevs))
                    else flags = combs[,byfacs] == bylevs
                    which(flags)
                })
                bylabs = apply(bygrid, 1, paste, collapse=",")
                bycols = sapply(byfacs, grep, names(combs))
                rnames = combs[ ,-bycols]
                if (!is.null(ncol(rnames))) rnames = apply(rnames, 1, paste, collapse=",")
            }
            
            # list to hold results
            Clist = list()
            zer = rep(0, ncol(K)) #### replaced nrow(lsms)) when lf arg added 
            
            # OK, let's go thru the bylist
            nby = length(bylist)
            for (i in 1:nby) {
                rows = bylist[[i]]
                cl = if(is.null(confcn)) contr[[method]] 
                    else confcn(rnames[rows] , ...)
                if (is.null(cl)) stop(paste("Unknown contrast family:", method))
                clx = lapply(cl, function(cc) {
                    ccc = zer; ccc[rows]=cc; ccc
                })
                if (nby > 1) names(clx) = paste(names(clx), "|", bylabs[i])
                
                # Currently, we're combining all these sets of contrasts in one table
                # Maybe we want to reconsider? If so, ship most of code below inside the
                # loop, and modify labels
                Clist = c(Clist, clx)
            }
            
            # Make a good label for the contrast table
            methdesc = attr(cl, "desc")
            
            # Try to at least explain why glht screws up
            if (!is.null(null.basis) && !is.null(glhargs)) {
                #glhargs = NULL
                warning("Error may occur in 'glht' due to rank deficiency")
            }
            if (lf || !is.null(glhargs)) { # create linear fcn for glht
                KK = t(sapply(Clist, function(con) {
                    nz = which(abs(con) > .0001)
                    K[ , nz] %*% con[nz]    
                }))
                if (lf) {
                    dimnames(KK)[[2]] = row.names(K)
                    ctbl = KK
                }
                else {
                    # If glht gets fixed for rank deficiency, may want to consider checking rows of KK
                    # for estimability (see code in do.est())
                    args = c(list(model=object, linfct=KK[ , used]), glhargs)
                    ctbl = summary(do.call("glht", args))
                }
            }
            else { # internal way of doing contrasts
                if (is.null(methdesc)) methdesc = method

                # Figure out the multiplicity adjustment
                adjattr = attr(cl, "adjust")
                if (autoadj) adj = ifelse(is.null(adjattr), no.adj, pmatch(adjattr, adjtbl))
                if (is.na(adj)) adj = no.adj
                
                ctbl = as.data.frame(t(sapply(Clist, function(con) {
                    nz = which(abs(con) > .0001)
                    k = K[ , nz] %*% con[nz]
                    do.est(k)
                })))
                
                # factors for mult adjustments...
                n.fam = nrow(lsms) / nby  ########### WAS sum(!is.na(lsms$lsmean)) / nby               
                n.contr = sum(!is.na(ctbl$estimate))
                if (!is.null(ddfm)) {
                    ctbl$t.ratio = round(ctbl$estimate / ctbl$SE, 5)
                    ctbl$p.value = round(adj.p.value(ctbl$t.ratio, ctbl$df, adj, n.fam, n.contr), 5)
                }
                else {
                    ctbl$z.ratio = round(ctbl$estimate / ctbl$SE, 5)
                    ctbl$p.value = round(adj.p.value(ctbl$z.ratio, 10000, adj, n.fam, n.contr), 5)
                }
                attr(ctbl, "mesg") = if(adj == 2)
                    paste("p values are adjusted using the", adjtbl[adj], "method for", n.fam, "means")
                else if (adj < length(adjtbl))
                    paste("p values are adjusted using the", adjtbl[adj], "method for", n.contr, "tests")
                else "p values are not adjusted"
                attr(ctbl, "print.row.names") = TRUE
                class(ctbl) = c("data.frame.lsm", "data.frame")
            }
            results[[paste(facs.lbl,methdesc)]] = ctbl
        }
    }
    if (!lf)
        class(results) = c("lsm","list")
    results
}


### S3 print method for "lsm" class - only reason we need this now is to support the 'omit' arg
print.lsm = function(x, omit=NULL, ...) {
    for (i in 1:length(x)) {
        if (i %in% omit) next
        cat(paste("$`", names(x)[i], "`\n", sep="")) # mimic print method for lists
        print(x[[i]])
        cat("\n")
    }
    invisible(x)
}

### S3 print method for annotated data.frames
print.data.frame.lsm = function(x, ...) {
    print.data.frame(x, row.names=attr(x, "print.row.names"))
    msg = attr(x, "mesg")
    if (!is.null(msg)) 
        for (j in 1:length(msg)) cat(paste("   ", msg[j], "\n"))
    invisible(x)
}

