##############################################################################
#    Copyright (c) 2012-2016 Russell V. Lenth                                #
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

# Namespace hooks for lsmeans

# .onLoad = function(libname, pkgname) {
#     # Function to check for existence of a variable
#     # This will be in the base package of R > 3.3.0
#     if (!exists("hasName", envir = getNamespace("utils"), inherits = FALSE)) {
#         assign("hasName", function(x, name)
#                 match(name, names(x), nomatch = 0L) > 0L,
#             envir = getNamespace(pkgname))
#     }
# }

# Just define the function for now. Maybe in a year we require R >= 3.4
hasName = function(x, name)
    match(name, names(x), nomatch = 0L) > 0L


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
        r = unique(ores[ores == paste0(flev, ".in.", nst)])
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

