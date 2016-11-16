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

### Support for overriding "ls" with "pm" in names

### general-purpose wrapper for creating pmxxxxx functions
.pmwrap = function(lsfcn, ...) {
    result = lsfcn(...)

        if (inherits(result, "ref.grid"))
        result = .sub.ls.pm(result)
    else if(inherits(result, "lsm.list")) {
        for (i in seq_along(result))
            result[[i]] = .sub.ls.pm(result[[i]])
        names(result) = gsub("^ls", "pm", names(result))
    }
    result
}

# returns an updated ref.grid or lsmobj with setName "ls..." replaced by "pm..."
.sub.ls.pm = function(object) {
    nm = object@misc$estName
    update(object, estName = gsub("^ls", "pm", nm))
}

### Exported implementations

pmmeans = function(...)
    .pmwrap(lsmeans, ...)

# uh, maybe not...   pmms = pmmeans

pmtrends = function(...)
    .pmwrap(lstrends, ...)


pmmip = function(...)
    lsmip(...)

pmm = function(...)
    lsm(...)

pmmobj = function(...)
    .pmwrap(lsmobj, ...)

pmm.options = function(...)
    lsm.options(...)

get.pmm.option = function(...)
    get.lsm.option(...)
