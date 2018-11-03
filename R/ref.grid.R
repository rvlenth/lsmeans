##############################################################################
#    Copyright (c) 2012-2018 Russell V. Lenth                                #
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

ref.grid = function(object, ...) 
    emmeans::ref_grid(object, ...)

recover.data = function(object, ...)
    emmeans::recover_data(object, ...)

lsm.basis = function(object, ...)
    emmeans::emm_basis(object, ...)


## methods for lsm.list

as.emm_list = function(object) {
    class(object) = c("emm_list", "list")
    object
}

as.glht.lsm.list  = function(object, ...) 
    emmeans::as.glht(as.emm_list(object), ...)

coef.lsm.list  = function(object, ...) 
    stats::coef(as.emm_list(object), ...)

confint.lsm.list  = function(object, ...) 
    stats::confint(as.emm_list(object), ...)

contrast.lsm.list  = function(object, ...) 
    emmeans::contrast(as.emm_list(object), ...)

pairs.lsm.list  = function(x, ...) 
    graphics::pairs(as.emm_list(x), ...)

summary.lsm.list = function(object, ...)
    summary(as.emm_list(object), ...)

str.lsm.list  = function(object, ...) 
    utils::str(as.emm_list(object), ...)

test.lsm.list  = function(object, ...) 
    emmeans::test(as.emm_list(object), ...)

