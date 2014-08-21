# S3 plot method for lsmobj objects (NOT ref.grid as relies on pri.fac attribute etc.)
# ... are arguments sent to update()


plot.lsmobj = function(x, y, type, comparisons = FALSE, 
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
    plot(summ, extra = extra, ...)
}

# May use in place of plot.lsmobj but no control over level etc.
# extra is a placeholder for comparison-interval stuff
plot.summary.ref.grid = function(x, y, horizontal = TRUE, xlab, ylab, extra = NULL, ...) {
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
    prepanel.ci = function(x, y, horizontal=TRUE, lcl, ucl, subscripts, ...) {
        x = as.numeric(x)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        if (horizontal)
            list(xlim = range(x, ucl, lcl, finite = TRUE)) 
        else
            list(ylim = range(y, ucl, lcl, finite = TRUE)) 
    }
    panel.ci <- function(x, y, horizontal=TRUE, lcl, ucl, locl, licl, ricl, rocl,
                         subscripts, pch = 16, 
                         lty = dot.line$lty, lwd = dot.line$lwd, 
                         col = dot.symbol$col, col.line = dot.line$col, ...) {
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        x = as.numeric(x)
        y = as.numeric(y)
        lcl = as.numeric(lcl[subscripts])
        ucl = as.numeric(ucl[subscripts])
        compare = !is.null(locl)
        if(compare) {
            locl = as.numeric(locl[subscripts])
            licl = as.numeric(licl[subscripts])
            ricl = as.numeric(ricl[subscripts])
            rocl = as.numeric(rocl[subscripts])
        }
        if(horizontal) {
            panel.abline(h = unique(y), col = col.line, lty = lty, lwd = lwd)
            panel.arrows(lcl, y, ucl, y, col = col, length = .6, unit = "char", angle = 90, code = 3)
            if(compare) {
                s = (x > min(x))
                panel.arrows(locl[s], y[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "gray", type = "closed", fill="gray")
                panel.arrows(licl[s], y[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                s = (x < max(x))
                panel.arrows(rocl[s], y[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "gray", type = "closed", fill="gray")
                panel.arrows(ricl[s], y[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "red", type = "closed", fill="red")
            }
        }
        else {
            panel.abline(v = unique(x), col = col.line, lty = lty, lwd = lwd)
            panel.arrows(x, lcl, x, ucl, col=col, length = .6, unit = "char", angle = 90, code = 3)
            if(compare) {
                s = (y > min(y))
                panel.arrows(x[s], locl[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "gray", type = "closed", fill="gray")
                panel.arrows(x[s], licl[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                s = (y < max(y))
                panel.arrows(x[s], rocl[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "gray", type = "closed", fill="gray")
                panel.arrows(x[s], ricl[s], x[s], y[s], length = .75, unit = "char", code = 1, col = "red", type = "closed", fill="red")
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
    
    # Obtain comparison limits - locl is left outer comparison limit, etc.
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
        involved = lapply(seq_len(neach), function(x) union(which(id2==x), which(id1==x)))
        mind = maxd = numeric(length(lbv))
        for (by in ubv) {
            d = del[pbv == by]
            rows = which(lbv == by)
            for(i in seq_len(neach)) {
                r = range(d[involved[[i]]])
                mind[rows[i]] = r[1]
                maxd[rows[i]] = r[2]
            }
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
        
        locl = invtran(est - maxd)
        licl = invtran(est - mind)
        ricl = invtran(est + mind)
        rocl = invtran(est + maxd)
    }
    else locl = licl = ricl = rocl = NULL
    
    
    
    facName = paste(priv, collapse=":")
    form = as.formula(chform)
    if (horizontal) {
        if (missing(xlab)) xlab = estName
        if (missing(ylab)) ylab = facName
        dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = TRUE,
                ylab = ylab, xlab = xlab,
                data = summ, lcl=lcl, ucl=ucl, 
                locl=locl, licl=licl, ricl=ricl, rocl=rocl, ...)
    }
    else {
        if (missing(xlab)) xlab = facName
        if (missing(ylab)) ylab = estName
        dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                strip = my.strip, horizontal = FALSE,
                xlab = paste(priv, collapse=":"), ylab = ylab,
                data = summ, lcl=lcl, ucl=ucl, 
                locl=locl, licl=licl, ricl=ricl, rocl=rocl, ...)
    }
    
    
}
