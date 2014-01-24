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

