lsmeans = function(object, specs, adjust=c("auto","tukey","sidak","scheffe",p.adjust.methods), conf = .95, 
                   at, trend, contr=list(), 
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
    
# added 6-18-2013 - allow cov.reduce to be logical
    if (is.logical(cov.reduce)) {
        if (cov.reduce) cov.reduce = function(x, name) mean(x)
        else cov.reduce = function(x, name) sort(unique(x))
    }
    
    # for later use
    adjtbl = c("auto","tukey","sidak","scheffe",p.adjust.methods)
    no.adj = pmatch("none", adjtbl)
    adj = pmatch(adjust, adjtbl)[1]
    if (is.na(adj)) {
        adj = 1
        warning("Unknown or non-unique `adjust' method -- automatic method will be used")
    }
    autoadj = (adj == 1)
    
    trend.flag = !missing(trend)
    
# Get model formula and index of response
    if (inherits(object, "gls")) {
        Terms = getCovariateFormula(object)
        yidx = 0
    }
    else {
        Terms = terms(object)
        yidx = attr(Terms, "response")
    }
### was    Terms = delete.response(terms(object))
### but we need to keep track of complete cases, including y values
### (Correction 2-13-13, version 1.06-05)    
    
    # get the pure formula w/o extra stuff
    formrhs = formula(Terms)
    
# ddfm will be replaced with a function of k and se if there is a way to get denom df    
    ddfm = adjV = NULL
    
# Figure out thecall (fixed effects part of model), bhat (=coefs), contrasts attr
    if (inherits(object, "mer") || inherits(object, "merMod")) {
        if(!isLMM(object) && !isGLMM(object)) 
            stop("Can't handle a nonlinear mixed model")
        thecall = slot(object, "call")
        bhat = fixef(object)
        contrasts = attr(model.matrix(object), "contrasts")
        if (isLMM(object)) {
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
        if (ncol(tR) < nrow(tR)) # add columns if not square
            tR = cbind(tR, matrix(0, nrow=nrow(tR), ncol=nrow(tR)-ncol(tR)))
        rank = object$qr$rank
        # last few rows are zero -- add a diagonal
        for (i in (rank+1):nrow(tR)) 
            tR[i,i] = 1
        null.basis = qr.resid(qr(tR[, seq_len(rank)]), tR[, -seq_len(rank)])
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
                    subset = eval(thecall$subset, enclos=envir),
                    na.action = na.omit, drop.unused.levels = TRUE)
### Correction 2-13-13, version 1.06-05 -- 
### added drop.unused terms and na.omit to above
    
    # Now X contains the data used to fit the model, w/o any expansions (e.g. poly() calls)

# Start accumulating info for the vars. 
# baselevs has the levels of all factors, or the "at" values for all covariates
# xlev has the factor levels only, for use in model.frame calls
    baselevs = xlev = matdat = list()
    # allow a vector of character strings
    if (is.character(specs)) specs = as.list(specs)
    # allow a single formula
    if (!is.list(specs)) specs = list(specs)
    
    all.var.names = names(X)
 
    for (xname in all.var.names) {
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
            # Keep track of covariates with more than one level
            #    if (length(baselevs[[xname]]) > 1)
            #        mult.covar = c(mult.covar, xname)
                
            }
        }
    }
    
    # If 'trend' present, we'll set up a difference quotient based on 
    # a fraction of the range of the variable
    if (trend.flag) {
        # trend could be a variable or a term; 
        # set up trend.xnm to always be variable names(s)
        if (! is.character(trend))
            stop("'trend' must be of character type")
        trend.xnm = trend[1]
        trend.x = X[[trend.xnm]]
        if (is.null(trend.x)) { 
# 'trend' not found, so hold out for possibility that it is the name of a term
            trend.xnm = try(all.vars(as.formula(paste("~",trend))), silent=TRUE)
            if (inherits(trend.xnm, "try-error")) 
                trend.xnm = character(0)
            trend.h = -1 # flag that trend is a term
        }
        else {
            if (!is.numeric(trend.x))
                stop("'trend' must refer to a numeric predictor")
            trend.h = diff(range(trend.x))*.001
            baselevs[[trend]] = baselevs[[trend]][1] + c(-1,1) * trend.h / 2
            # rearrange ordering so trend variable is first - makes bookkeeping easier
            sidx = match(trend, names(baselevs))
            baselevs = c(baselevs[sidx], baselevs[-sidx])
        }
    }
    
    # Keep the response variable from enlarging the grid, no matter what
    if (yidx > 0) {
        yname = as.character(attr(Terms, "variables")[[1 + yidx]])
        if (!is.na(match(yname, names(baselevs))[1])) 
            baselevs[[yname]] = NA
    }
    
    # OK. Now make a grid of the factor levels of interest, along w/ covariate "at" values
    grid = do.call(expand.grid, baselevs)
    
    # add any matrices
    for (nm in names(matdat))
        grid[[nm]] = matrix(rep(matdat[[nm]], each=nrow(grid)), nrow=nrow(grid))

    # It turns out that numerics coerced to factors are a real pain in the butt when it comes
    # to matching levels. Life will be simpler if we turn them into factors in the X matrix 
    # and update the base levels accordingly with the same labels
    #
    #--- Version 1.10 - I don't think I need this anymore with new matching routine
#     for (var in coerced) {
#         X[[var]] = factor(X[[var]])
#         baselevs[[var]] = levels(X[[var]])
#     }

    # Now make a new dataset with just the factor combs and covariate values we want for prediction
    # WARNING -- This will overwrite X, so get anything you need from X BEFORE we get here
    m = model.frame(Terms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(Terms, m, contrasts.arg = contrasts)

    # Compute the derivatives when 'trend' is provided
    if(trend.flag) {
        if (trend.h < 0) { # 'trend' is a term - fix the coefs
            # wipe-out spaces to make it less fussy
            term.nm = gsub(" ", "", dimnames(X)[[2]])
            trend = gsub(" ", "", trend)
            trend.idx = match(trend, term.nm)
            if (is.na(trend.idx))
                stop(paste("trend value '", trend, "' is neither a variable nor a model term", sep=""))
            prev.X = X
# Symbolic differentiation... Replace columns of X with:
#  1 if it is col for trend
#  0 if it does not contain trend
#  previous version of X[,j] where j is interaction of other predictors            
            for (i in seq_len(ncol(X))) {
                trm = term.nm[i]
                trm.pieces = strsplit(trm, ":")[[1]]
                trm.mat = match(trend, trm.pieces)
                if (is.na(trm.mat))
                    X[, i] = 0
                else {
                    trm.otr = paste(trm.pieces[-trm.mat], collapse=":")
                    trm.ref = match(trm.otr, term.nm)
                    if (is.na(trm.ref))     X[, i] = 0 + (i == trend.idx)
                    else                    X[, i] = prev.X[, trm.ref]
                }
            }
        }
        else { # 'trend' is a variable - do the diff quotient
            evens = 2 * (seq_len(nrow(X)/2))
            X = (X[evens, ] - X[evens-1, ]) / trend.h
            # restore the base levels
            baselevs[[trend]] = baselevs[[trend]][1] + trend.h/2
            grid = grid[evens, , drop=FALSE]
        }
    }
    
    # If necessary revise grid with coerced numeric factors replaced with factor levels
    if (length(coerced) > 0) grid = do.call("expand.grid", baselevs)
    
    # All factors (excluding covariates)
    # version 1.10 - no longer excluding covariates
    allFacs = all.var.names
    
    ### Array of indexes for rows of X, organized by dimensions
    row.indexes = array(seq_len(nrow(X)), sapply(baselevs, length))
    
    
    
    # Get a vector of terms in the model, for checking
    mod.terms = attr(Terms, "term.labels")
        ### was strsplit(as.character(formrhs[2])[[1]], "\\+")[[1]]
    
##### routine returns TRUE iff all elements of facs are contained in a model term with another predictor
    some.term.contains = function(facs) {
        # When checking trend, add that variable to the mix
        if(trend.flag) 
            facs = union(facs, trend.xnm)
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
    for (i in seq_len(length(specs))) {
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

### New (version 1.10) more efficient derivation of K matrix
        RI = plyr:::splitter_a(row.indexes, match(facs, names(baselevs)))
    # Each entry of RI has the row indexes of X
    # for each combination of facs (in expand.grid order)
        K = sapply(RI, function(idx) {
            fac.reduce(X[idx, , drop=FALSE], "")
        })
                
#--- above code replaces pre-1.10 code below...
#         # For each comb, find the needed lin. comb. of bhat to estimate
#         # (These will end up being the COLUMNS of K)
#         K = apply(combs, 1, function(lev) {
#             matches = apply(grid, 1, function(row) {
#                 if (is.numeric(lev)) 
#                     all(zapsmall(as.numeric(row[facs]) - lev) == 0) 
#                 else
#                     all(row[facs] == lev)
#             })
#             nmat = sum(matches)
#             if (nmat == 0) stop(paste("Can't predict at level", lev, "of", "facs.lbl"))
#             else fac.reduce(X[matches, , drop=FALSE], lev)
#         })
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
    # I added zapsmall in ptukey call, seems to help with its flaky behavior
        adj.p.value = function(t, df, meth, fam.size, n.contr) {
            abst = abs(t)
            if (meth <= 4)
                switch(meth,
                       NA,                                                     # should not happen
                       ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),   # tukey
                       1 - (1 - 2*pt(abst, df, lower.tail=FALSE))^n.contr,     # sidak
                       pf(t^2/(fam.size-1), fam.size-1, df, lower.tail=FALSE)  # scheffe
                )
            else
                p.adjust(2*pt(abst, df, lower.tail=FALSE), adjtbl[meth], n=n.contr)
        }
        
        # LS means
        if (!trend.flag) {
            effname = "lsmean"
            lsmentry = paste(facs.lbl, "lsmeans")
        }
        else {
            effname = paste(trend, "trend", sep=".")
            lsmentry = paste(effname, "by", facs.lbl)
        }

        if (lf) {
            results[[lsmentry]] = t(K)
        }
        else {
            lsms = as.data.frame(t(apply(K,2,do.est)))
            # fix-up names and get CIs
            names(lsms)[1] = effname
            # include factor levels
            lsms = cbind(combs, lsms)
            if (conf > 1) conf = conf/100 # pct --> frac
            if ((conf < 1) && (conf > .01)) {
                if (is.null(ddfm)) {
                    me = qnorm((1-conf)/2, lower.tail=FALSE) * lsms$SE
                    lsms$asymp.LCL = lsms[[effname]] - me
                    lsms$asymp.UCL = lsms[[effname]] + me
                }
                else {
                    me = qt((1-conf)/2, lsms$df, lower.tail=FALSE) * lsms$SE
                    lsms$lower.CL = lsms[[effname]] - me
                    lsms$upper.CL = lsms[[effname]] + me
                }
            }
            attr(lsms, "print.row.names") = FALSE
            class(lsms) = c("data.frame.lsm", "data.frame")
            results[[lsmentry]] = lsms
        }        
        
        # Do requested contrasts
        if (! is.null(method)) {
            # look for contrast fcn
            fn = paste(method, "lsmc", sep=".")
            confcn = if (exists(fn, mode="function")) get(fn) 
                else NULL
            
            # bylist will be a list of subsets of the combs to be contrasted
            if (is.null(byfacs)) bylist = list(seq_len(nrow(combs))) # all in one set
            else {
                bg = list()
                for (f in byfacs) bg[[f]] = baselevs[[f]]
                bygrid = do.call("expand.grid", bg)
                bylist = lapply(seq_len(nrow(bygrid)), function(row) {
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
            for (i in seq_len(nby)) {
                rows = bylist[[i]]
                cl = if(is.null(confcn)) contr[[method]] 
                    else confcn(rnames[rows] , ...)
                if (is.null(cl)) stop(paste("Unknown contrast family:", method))
                clx = lapply(cl, function(cc) {
                    if(length(cc) != length(rows))
                        stop(paste(length(cc), " contrast coefficients in '", method, 
                                   "' when ", length(rows), " were expected", sep=""))
                    ccc = zer; 
                    ccc[rows]=cc; 
                    ccc
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

#---removed - not really needed, & apparently not more efficient than splitter_a
# # My version of plyr:::splitter_a, but returns a matrix
# # On return, each column has the elements of .array for each combination of .margins
# # Order of columns is same as obtained using expand.grid with the same variables
# .mysplit = function(.array, .margins) {
#     dims = dim(.array)
#     len = length(dims)
#     if (any(.margins > len))
#         stop ("'.margins' exceeds dimensions of '.array'")
#     prm = c(setdiff(seq_len(len), .margins), .margins)
#     matrix(aperm.default(.array, prm, FALSE), ncol = prod(dims[.margins]))
# }

### S3 print method for "lsm" class - only reason we need this now is to support the 'omit' arg
print.lsm = function(x, omit=NULL, ...) {
    for (i in seq_len(length(x))) {
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
        for (j in seq_len(length(msg))) cat(paste("   ", msg[j], "\n"))
    invisible(x)
}

