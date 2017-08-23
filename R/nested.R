##############################################################################
#    Copyright (c) 2012-2017 Russell V. Lenth                                #
#                                                                            #
#    This file is part of the lsmeans package for R (*lsmeans*)              #
#                                                                            #
#    *lsmeans* is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 2 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    *lsmeans* is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with R and *lsmeans*.  If not, see                                #
#    <https://www.r-project.org/Licenses/> and/or                            #
#    <http://www.gnu.org/licenses/>.                                         #
##############################################################################

# Code supporting nested models

# This code replies on nested structures specified in a named list like
#    list(a = "b", c = c("d", "e"))
# ... to denote a %in% b, c %in% d*e

# Internal function to deal with nested structures. 
#   rgobj        -- a ref.grid
#   specs, ...   -- arguments for lsmeans
#   nesting      -- a named list of nesting info
# This function works by subsetting rgobj as needed, and applying lsmeans
# to each subsetted object
# This is a servant to lsmeans.character.ref.grid, so we can assume specs is character
.nested_lsm = function(rgobj, specs, by = NULL, ..., nesting) {
    # # Trap something not supported for these... This doesn't work
    # dots = list(...)
    # if("weights" %in% dots)
    #     if(!is.na(pmatch(dots$weights, "show.levels")))
    #         stop('weights = "show.levels" is not supported for nested models.')

    #### Two issues to worry about....
    # (1) specs contains nested factors. We need to include their grouping factors
    xspecs = intersect(union(specs, by), names(nesting))
    if (length(xspecs) > 0) {
        xgrps = unlist(nesting[xspecs])
        specs = union(union(xspecs, xgrps), specs)  # expanded specs with flagged ones first
        by = setdiff(by, xspecs) # can't use nested factors for grouping
    }
    # (2) If we average over any nested factors, we need to do it separately
    avg.over = setdiff(names(rgobj@levels), union(specs, by))
    afacs = intersect(names(nesting), avg.over) ### DUH!names(nesting)[names(nesting) %in% avg.over]
    if (length(afacs) == 0)  { # no nesting issues; just use lsmeans
        result = lsmeans(rgobj, specs, by = by, ...)
    }
    else { # we need to handle each group separately
        sz = sapply(afacs, function(nm) length(nesting[[nm]]))
        # use highest-order one first: potentially, we end up calling this recursively
        afac = afacs[rev(order(sz))][1] 
        grpfacs = nesting[[afac]]
        gspecs = union(specs, union(by, grpfacs))
        grpids = as.character(interaction(rgobj@grid[, grpfacs]))
        grps = do.call(expand.grid, rgobj@levels[grpfacs])  # all combinations of group factors
        result = NULL
        rg = rgobj
        for (i in seq_len(nrow(grps))) {
            sig = as.character(interaction(grps[i, ]))
            rows = which(grpids == sig)
            grd = rgobj@grid[rows, , drop = FALSE]
            lf = rgobj@linfct[rows, , drop = FALSE]
            # Reduce grid to infacs that actually appear in  this group
            nzg = grd[grd$.wgt. > 0, , drop = FALSE]
            rows = integer(0)
            # focus on levels of afac that exist in this group
            levs = unique(nzg[[afac]])
            rg@levels[[afac]] = levs
            rows = union(rows, which(grd[[afac]] %in% levs))
            rg@grid = grd[rows, , drop = FALSE]
            rg@linfct = lf[rows, , drop = FALSE]
            for (j in seq_along(grpfacs))
                rg@levels[[grpfacs[j]]] = grps[i, j]
            lsm = suppressMessages(lsmeans(rg, gspecs, ...))
            if (is.null(result))
                result = lsm
            else {
                result@grid = rbind(result@grid, lsm@grid)
                result@linfct = rbind(result@linfct, lsm@linfct)
            }
        }
        for (j in seq_along(grpfacs))
            result@levels[grpfacs[j]] = rgobj@levels[grpfacs[j]]
        
        result@misc$avgd.over = setdiff(union(result@misc$avgd.over, avg.over), gspecs)
        nkeep = intersect(names(nesting), names(result@levels))
        if (length(nkeep) > 0)
            result@model.info$nesting = nesting[nkeep]
        
        # Note: if any nesting remains, this next call recurs back to this function
        result = lsmeans(result, specs, by = by, ...)
    }
    
    if (length(xspecs) > 0)
        result@misc$display = .find.nonempty.nests(result, xspecs, nesting)
    else
        result@misc$display = NULL
    
    # preserve any nesting that still exists
    result@model.info$nesting = nesting[names(nesting) %in% names(result@levels)]
    result
}


### contrast function for nested structures
.nested_contrast = function(rgobj, method = "eff", by = NULL, adjust, ...) {
    nesting = rgobj@model.info$nesting
    # Prevent meaningless cases -- if A %in% B, we can't have A in 'by' without B
    # Our remedy will be to EXPAND the by list
    for (nm in intersect(by, names(nesting)))
        if (!all(nesting[[nm]] %in% by)) {
            by = union(by, nesting[[nm]])
            message("Note: Grouping factor(s) for '", nm, "' have been added to the 'by' list.")
        }
    facs = setdiff(names(nesting), by)
    if (length(facs) == 0)
        stop("There are no factor levels left to contrast. Try taking nested factors out of 'by'.")
    
    if(!is.character(method))
        stop ("Non-character contrast methods are not supported with nested objects")
    
    testcon = get(paste0(method, ".lsmc"))(1:3)
    if(missing(adjust)) 
        adjust = attr(testcon, "adjust")
    estType = attr(testcon, "type")

    wkrg = rgobj # working copy
    facs = setdiff(names(wkrg@levels), by)  # these are the factors we'll combine & contrast
    if (!is.null(display <- wkrg@misc$display))
        wkrg = wkrg[which(display), drop.levels = TRUE]
    wkrg@model.info$nesting = wkrg@misc$display = NULL
    by.rows = .find.by.rows(wkrg@grid, by)
    if(length(by.rows) == 1)
        result = contrast(wkrg, method = method, by = by, ...)
    else {
        result = lapply(by.rows, function(rows) {
            contrast.ref.grid(wkrg[rows, drop.levels = TRUE], method = method, 
                              by = by, adjust = adjust, ...)
        })
        result$adjust = ifelse(is.null(adjust), "none", adjust)
        result = do.call(rbind.ref.grid, result)
        result = update(result, by = by, 
                        estType = ifelse(is.null(estType), "contrast", estType))
        cname = setdiff(names(result@levels), by)
        result@model.info$nesting[[cname]] = by
    }
    result@misc$orig.grid = result@misc$con.code = NULL

    for (nm in by) {
        if (nm %in% names(nesting))
            result@model.info$nesting[[nm]] = intersect(nesting[[nm]], by)
    }
    result
}


# Internal function to find nonempty cells in nested structures in rgobj for xfacs
# Returns logical vector, FALSE are rows of the grid we needn't display
.find.nonempty.nests = function(rgobj, xfacs, nesting = rgobj@model.info$nesting) {
    grid = rgobj@grid
    keep = rep(TRUE, nrow(grid))
    for (x in xfacs) {
        facs = union(x, nesting[[x]])
        combs = do.call(expand.grid, rgobj@levels[facs])
        levs = as.character(interaction(combs))
        glevs = as.character(interaction(grid[facs]))
        
        for (lev in levs) {
            idx = which(glevs == lev)
            if (all(grid$.wgt.[idx] == 0)) {
                keep[idx] = FALSE
                levs[levs==lev] = ""
            }
        }
    }
    keep
}


# Internal function to find nesting
# We look at two things:
# (1) structural nesting - i.e., any combinations of
#     factors A and B for which each level of A occurs with one and only one
#     level of B. If so, we deem A %in% B.
# (2) Model-term nesting - cases where a factor appears not as a main effect
#     but only in higher-order terms. This is discovered using the 1s and 2s in 
#     trms$factors
# The function returns a named list, e.g., list(A = "B")
# If none found, an empty list is returned.
.find_nests = function(grid, trms) {
    result = list()
    nms = names(grid)[sapply(grid, is.factor)]
    if (length(nms) < 2)
        return (result)
    g = grid[grid$.wgt. > 0, nms, drop = FALSE]
    for (nm in nms) {
        x = levels(g[[nm]])
###        x = x[x != ".nref."]     # ignore reference level from nested()
        otrs = nms[!(nms == nm)]
        max.levs = sapply(otrs, function(n) {
            max(sapply(x, function(lev) length(unique(g[[n]][g[[nm]] == lev]))))
        })
        if (any(max.levs == 1))
            result[[nm]] = otrs[max.levs == 1]
    }
    
    # Now look at factors attribute
    fac = attr(trms, "factors")
    fac = fac[intersect(nms, row.names(fac)), , drop = FALSE]
    for (j in seq_len(ncol(fac))) {
        if (any(fac[, j] == 2)) {
            nst = nms[fac[, j] == 1]
            for (nm in nst)
                result[[nm]] = nms[fac[, j] == 2]
        }
    }
    
    result
}

# internal function to format a list of nested levels
.fmt.nest = function(nlist) {
    if (length(nlist) == 0)
        "none"
    else {
        tmp = lapply(nlist, function(x) paste(x, collapse = "*"))
        paste(sapply(names(nlist), function (nm) paste0(nm, " %in% ", tmp[[nm]])),
              collapse = ", ")
    }
}


# ### I'm removing this because I now think it creates more problems than it solves
# #
# # courtesy function to create levels for a nested structure factor %in% nest
# # factor: factor (or interaction() result)
# # ...:    factors in nest
# # SAS:    if (FALSE|TRUE), reference level in each nest is (first|last)
# nested = function(factor, ..., SAS = FALSE) {
#     nfacs = list(...)
#     if (length(nfacs) == 0)
#         return(factor)
#     nfacs$drop = TRUE
#     nest = do.call(interaction, nfacs)
#     result = as.character(interaction(factor, nest, sep = ".in."))
#     ores = unique(sort(result))
#     nlev = levels(nest)
#     flev = levels(factor)
#     refs = lapply(nlev, function(nst) {
#         r = ores[ores %in% paste0(flev, ".in.", nst)]
#         ifelse (SAS, rev(r)[1], r[1])
#     })
#     result[result %in% refs] = ".nref."
#     ores[ores %in% refs] = ".nref."
#     ores = setdiff(ores, ".nref.")
#     if (SAS)
#         factor(result, levels = c(ores, ".nref."))
#     else
#         factor(result, levels = c(".nref.", ores))
# }

