# Namespace hooks for lsmeans

.onLoad = function(libname, pkgname) {
    # Function to check for existence of a variable
    # This will be in the base package of R > 3.3.0
    if (!exists("hasName", envir = getNamespace("utils"), inherits = FALSE)) {
        assign("hasName", function(x, name) 
                match(names(x), name, nomatch = 0L) > 0L,
            envir = getNamespace(pkgname))
    }
}
