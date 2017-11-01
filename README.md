R package **lsmeans**: Least-squares means (estimated marginal means)
====

[![cran version](http://www.r-pkg.org/badges/version/lsmeans)](https://cran.r-project.org/package=lsmeans)
[![downloads](http://cranlogs.r-pkg.org/badges/lsmeans)](http://cranlogs.r-pkg.org/badges/lsmeans)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/lsmeans)](http://cranlogs.r-pkg.org/badges/grand-total/lsmeans)
[![Research software impact](http://depsy.org/api/package/cran/lsmeans/badge.svg)](http://depsy.org/package/r/lsmeans)

## NEWS
The **lsmeans** package is being deprecated. Users are encouraged to
switch to [**emmeans**](https://github.com/rvlenth/emmeans)
(estimated marginal means), now available on CRAN.
The **lsmeans** package will be archived on CRAN at some not-too-distant
time in the future. Note that:

  * R scripts that use **lsmeans** will still work with **emmeans** after making 
    minor changes (use `emmeans:::convert_scripts()`). 
  * Existing objects created with **lsmeans** can be converted to work 
    with the new package via `emmeans:::convert_workspace()`. 
  * See `vignette("transition-from-lsmeans", "emmeans")` for more details.
  