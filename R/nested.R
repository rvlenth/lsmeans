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


# Internal function to deal with nested structures. 
#   rgobj        -- a ref.grid
#   specs, ...   -- arguments for lsmeans
#   nests        -- a formula or list of formulas each of the form  factors ~ in.factors
# This function works by subsetting rgobj as needed, and applying lsmeans
# to each subsetted object
# This is a servant to lsmeans.character.ref.grid, so we can assume specs is character
.nested_lsm = function(rgobj, specs, by = NULL, ..., nests) {
    if (inherits(nests, "formula"))   # for a single nested structure
        nests = list(nests)
    avg.over = setdiff(names(rgobj@levels), union(specs, by))
    infacs = grpfacs = character(0)
    # figure out which nested factors get averaged over, and which groups are involved
    for (nest in nests) { 
        if (!inherits(nest, "formula") || length(nest) < 3)
            stop("All nesting specifications must be two-sided formulas")
        gf = .all.vars(nest[-2])
        nf = .all.vars(nest[-3])
        af = nf[nf %in% avg.over]
        if (length(af) > 0) {
            infacs = union(infacs, af)
            grpfacs = union(grpfacs, gf)
        }
    }
    if (length(grpfacs) == 0)  # no nesting issues; just use lsmeans
        result = lsmeans(rgobj, specs, ...)
    else { # we need to handle each group separately
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
            for (fac in infacs) {
                levs = unique(nzg[[fac]])
                rg@levels[[fac]] = levs
                rows = union(rows, which(grd[[fac]] %in% levs))
            }
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
        # # Take care of a subtlety: What we have constructed always will have groups 
        # # as the last sets of levels
        # lvls = names(result@levels)
        # newlvls = union(setdiff(lvls, grpfacs), grpfacs)
        # result@grid = result@grid[, newlvls, drop = FALSE]
        # result@levels = result@levels[newlvls]
        result@misc$avgd.over = setdiff(union(result@misc$avgd.over, avg.over), gspecs)
        result = lsmeans(result, specs, by = by, ...)
    }
    result
}



# courtesy function to create levels for a nested structure factor %in% nest
# factor: factor (or interaction() result)
# ...:    factors in nest
# SAS:    if (FALSE|TRUE), reference level in each nest is (first|last)
nested = function(factor, ..., SAS = FALSE) {
    nfacs = list(...)
    if (length(nfacs) == 0)
        return(factor)
    nfacs$drop = TRUE
    nest = do.call(interaction, nfacs)
    result = as.character(interaction(factor, nest, sep = ".in."))
    ores = unique(sort(result))
    nlev = levels(nest)
    flev = levels(factor)
    refs = lapply(nlev, function(nst) {
        r = ores[ores %in% paste0(flev, ".in.", nst)]
        ifelse (SAS, rev(r)[1], r[1])
    })
    result[result %in% refs] = "ref"
    ores[ores %in% refs] = "ref"
    ores = setdiff(ores, "ref")
    if (SAS)
        factor(result, levels = c(ores, "ref"))
    else
        factor(result, levels = c("ref", ores))
}

