# Namespace hooks for lsmeans

.onLoad = function(libname, pkgname) {
    # Function to check for existence of a variable
    # This will be in the base package of R > 3.3.0
    if (!exists("hasName", envir = getNamespace("base"), inherits = FALSE)) {
        assign("hasName", function(x, name) name %in% names(x),
               envir = getNamespace(pkgname))
    }
}
