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

# interaction-plot code
#
# This is just sketched out, but the idea is that we have generic
# int.plot, for which the S3 method for 'object' of class data.frame 
# becomes the generic for the second argument 'fac'

int.plot = function(object, ...)
    UseMethod("int.plot")

# S3 method for data.frame becomes generic for second arg
int.plot.data.frame = function(object, fac, ...)
    UseMethod("int.plot.data.frame", fac)


# For a model object... use lsmeans
int.plot.default = function(object, fac, ...) {
    df = summary(lsmeans(object, ~Variety)) # need to extract OK arguments though
    int.plot(df, fac, ...)
}
    
int.plot.data.frame.default = function(object, fac, ...) {
    cat("Reached int.plot.data.frame.default\n")
}
int.plot.data.frame.formula = function(object, fac, ...) {
    cat("Reached int.plot.data.frame.formula\n")
}
int.plot.data.frame.character = function(object, fac, ...) {
    cat("Reached int.plot.data.frame.character\n")
}

