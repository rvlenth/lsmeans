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

# Methods for lsm.list objects

# Summary method for an lsm.list
summary.lsm.list <- function(object, ...)
    lapply(object, function(x) {
        if (inherits(x, "summary.ref.grid"))  x
        else summary(x, ...)
    })

print.lsm.list <- function(x, ...) 
    print(summary(x, ...))

str.lsm.list = function(object, ...) {
    for(nm in names(object)) {
        cat(paste("$", nm, "\n", sep=""))
        str(object[[nm]])
        cat("\n")
    }
}

# Courtesy methods to make it more friendly for follow-ups
contrast.lsm.list = function(object, ... , which = 1) {
    contrast(object[[which]], ...)
}

pairs.lsm.list = function(x, ..., which = 1) {
    pairs(x[[which]], ...)
}

test.lsm.list = function(object, ..., which = 1) {
    test(object[[which]], ...)
}

confint.lsm.list = function(object, ..., which = 1) {
    confint(object[[which]], ...)
}

cld.lsm.list = function(object, ..., which = 1) {
    cld(object[[which]], ...)
}


