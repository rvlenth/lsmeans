# S3 plot method for lsmobj objects (NOT ref.grid as relies on pri.fac attribute etc.)
# ... are arguments sent to update()


plot.lsmobj = function(x, y, type, intervals = TRUE, comparisons = FALSE, 
                       c.alpha = .05, c.adjust = "tukey", ...) {
    if(!missing(type))
        object = update(x, predict.type = type, ..., silent = TRUE)
    else
        object = update(x, ..., silent = TRUE)
    summ = summary(object, infer = c(TRUE, FALSE))
    estName = attr(summ, "estName")
    extra = NULL
    if(comparisons) {
        extra = object
        extra@misc$comp.alpha = c.alpha
        extra@misc$comp.adjust = c.adjust
    }
    .plot.srg(x=summ, intervals = intervals, extra = extra, ...)
}

# May use in place of plot.lsmobj but no control over level etc.
# extra is a placeholder for comparison-interval stuff
plot.summary.ref.grid = function(x, y, horizontal = TRUE, xlab, ylab, ...) {
    .plot.srg (x, y, horizontal, xlab, ylab, ...)
}

# Workhorse for plot.summary.ref.grid
.plot.srg = function(x, y, horizontal = TRUE, xlab, ylab, intervals = TRUE, extra = NULL, ...) {
        
    if (!require("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    
    summ = x # so I don't get confused
    estName = attr(summ, "estName")
    clNames = attr(summ, "clNames")
    if (is.null(clNames)) {
        warning("No information available to display confidence limits")
        lcl = ucl = summ[[estName]]
    }
    else {
        lcl = summ[[clNames[1]]]
        ucl = summ[[clNames[2]]]
    }
    
    # Panel functions...
    prepanel.ci = function(x, y, horizontal=TRUE, intervals=TRUE,
                           lcl, ucl, subscripts, ...) {
        x = as.numeric(x)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        if (!intervals) # no special scaling needed
            list()
        else if (horizontal)
            list(xlim = range(x, ucl, lcl, finite = TRUE)) 
        else
            list(ylim = range(y, ucl, lcl, finite = TRUE)) 
    }
    panel.ci <- function(x, y, horizontal=TRUE, intervals=TRUE,
                         lcl, ucl, lcmpl, rcmpl,                          subscripts, pch = 16, 
                         lty = dot.line$lty, lwd = dot.line$lwd, 
                         col = dot.symbol$col, col.line = dot.line$col, ...) {
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        x = as.numeric(x)
        y = as.numeric(y)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        compare = !is.null(lcmpl)
        if(compare) {
            lcmpl = as.numeric(lcmpl[subscripts])
            rcmpl = as.numeric(rcmpl[subscripts])
        }
        if(horizontal) {
            panel.abline(h = unique(y), col = col.line, lty = lty, lwd = lwd)
            if(intervals) 
                panel.arrows(lcl, y, ucl, y, col = col, length = .6, unit = "char", angle = 90, code = 3)
            if(compare) {
                s = (x > min(x))
                panel.arrows(lcmpl[s], y[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                s = (x < max(x))
                panel.arrows(rcmpl[s], y[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
            }
        }
        else {
            panel.abline(v = unique(x), col = col.line, lty = lty, lwd = lwd)
            if(intervals)
                panel.arrows(x, lcl, x, ucl, col=col, length = .6, unit = "char", angle = 90, code = 3)
            if(compare) {
                s = (y > min(y))
                panel.arrows(x[s], lcmpl[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                s = (y < max(y))
                panel.arrows(x[s], rcmpl[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
            }
        }
        panel.xyplot(x, y, pch=16, ...)
    }
    my.strip = strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
    
    priv = attr(summ, "pri.vars")
    pf = do.call(paste, summ[priv])
    summ$pri.fac = factor(pf, levels=unique(pf))
    chform = ifelse(horizontal,
                    paste("pri.fac ~", estName),
                    paste(estName, "~ pri.fac"))
    
    byv = attr(summ, "by.vars")
    if (!is.null(byv)) {
        chform = paste(chform, "|", paste(byv, collapse="*"))
    }
    
    # Obtain comparison limits
    if (!is.null(extra)) {
        # we need to work on the linear predictor scale
        # typeid = 1 -> response, 2 -> other
        typeid = pmatch(extra@misc$predict.type, "response", nomatch = 2)
        if(length(typeid) < 1) typeid = 2        
        if (typeid == 1)
            est = predict(extra, type = "lp")
        else
            est = summ[[estName]]
        
        alpha = extra@misc$comp.alpha
        adjust = extra@misc$comp.adjust
        psumm = confint(pairs(extra), level = 1 - alpha, type = "lp", adjust = adjust)
        k = ncol(psumm)
        del = (psumm[[k]] - psumm[[k-1]]) / 4 # half the halfwidth, on lp scale
        diff = psumm[[attr(psumm, "estName")]]
        overlap = apply(psumm[ ,(k-1):k], 1, function(x) 2*min(-x[1],x[2])/(x[2]-x[1]))
        
        # figure out by variables and indexes
        if(is.null(byv)) {
            lbv = rep(1, nrow(summ))
            pbv = rep(1, nrow(psumm))
            ubv = 1
        }
        else {
            lbv = do.call("paste", summ[byv]) # strings for matching by variables
            pbv = do.call("paste", psumm[byv])
            ubv = unique(lbv)
        }
        neach = length(lbv) / length(ubv)
        # indexes for pairs results -- est[id1] - est[id2]
        id1 = rep(seq_len(neach-1), rev(seq_len(neach-1)))
        id2 = unlist(sapply(seq_len(neach-1), function(x) x + seq_len(neach-x)))
        # list of psumm row numbers involved in each summ row
        involved = lapply(seq_len(neach), function(x) union(which(id2==x), which(id1==x)))
        
        # initialize arrays
        mind = numeric(length(lbv))   # for minima of del
        llen = rlen = numeric(neach)  # for left and right arrow lengths
        npairs = length(id1)
        iden = diag(rep(1, 2*neach))
        
        for (by in ubv) {
            d = del[pbv == by]
            rows = which(lbv == by)
            for(i in seq_len(neach)) 
                mind[rows[i]] = min(d[involved[[i]]])
            
            # Set up regression equations to match arrow overlaps with interval overlaps
            # We'll add rows later (with weights 1) to match with mind values
            lmat = rmat = matrix(0, nrow = npairs, ncol = neach)
            y = numeric(npairs)
            v1 = 1 - overlap[pbv == by]
            dif = diff[pbv == by]
            for (i in 1:npairs) {
                wgt = 4 * max(0, ifelse(v1[i] < 1, v1[i], 2-v1[i]))
                # really this is sqrt of weight
                if (dif[i] > 0)   # id2  <----->  id1
                    lmat[i, id1[i]] = rmat[i, id2[i]] = wgt*v1[i]
                else  # id1  <----->  id2
                    rmat[i, id1[i]] = lmat[i, id2[i]] = wgt*v1[i]
                y[i] = wgt * abs(dif[i])
            }
            X = rbind(cbind(lmat, rmat),iden)
            y = c(y, rep(mind[rows], 2))
            soln = qr.coef(qr(X), y)
            llen[rows] = soln[seq_len(neach)]
            rlen[rows] = soln[neach + seq_len(neach)]
            
            # Perhaps put some kind of a check here?
        }
        invtran = I
        if (typeid == 1) {
            tran = extra@misc$tran
            if(is.character(tran)) {
                link = try(make.link(tran), silent=TRUE)
                if (!inherits(link, "try-error"))
                    invtran = link$linkinv
            }
            else if (is.list(tran))
                invtran = tran$linkinv
        }
        
        lcmpl = invtran(est - llen)
        rcmpl = invtran(est + rlen)
    }
    else lcmpl = rcmpl = NULL
    
    
    
    facName = paste(priv, collapse=":")
    form = as.formula(chform)
    if (horizontal) {
        if (missing(xlab)) xlab = estName
        if (missing(ylab)) ylab = facName
        dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = TRUE,
                ylab = ylab, xlab = xlab,
                data = summ, intervals = intervals, lcl=lcl, ucl=ucl, 
                lcmpl=lcmpl, rcmpl=rcmpl, ...)
    }
    else {
        if (missing(xlab)) xlab = facName
        if (missing(ylab)) ylab = estName
        dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = FALSE,
                xlab = paste(priv, collapse=":"), ylab = ylab,
                data = summ, intervals = intervals, lcl=lcl, ucl=ucl, 
                lcmpl=lcmpl, rcmpl=rcmpl, ...)
    }
}
