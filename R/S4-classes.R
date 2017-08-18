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

### S4 class definitions for lsmeans package


### ref.grid object -- for a reference grid
setClass("ref.grid", slots = c(
    model.info = "list",
    roles = "list",
    grid = "data.frame", 
    levels = "list",
    matlevs = "list",
    linfct = "matrix",
    bhat = "numeric",
    nbasis = "matrix",
    V = "matrix",
    dffun = "function",
    dfargs = "list",
    misc = "list",
    post.beta = "matrix"
))
# Note: misc will hold various extra params,
# including at least the following req'd by the summary method
#   estName: column name for the estimate in the summary ["prediction"]
#   infer: booleans (CIs?, tests?)  [(FALSE,FALSE)]
#   level: default conf level [.95]
#   adjust: default adjust method ["none"]
#   famSize: number of means in family



### lsmobj class -- almost trivial ext of ref.grid, structurally
# But origin can be very different from those of a reference grid
# In general its 'grid' will correspond to some set of 
# linear functions of grid points
setClass("lsmobj", contains="ref.grid")


