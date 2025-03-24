R package **lsmeans**: Least-squares means (estimated marginal means)
====

## NEWS
The **lsmeans** package is now just a front-end to **emmeans**. Users are encouraged to
switch to [**emmeans**](https://github.com/rvlenth/emmeans)
(estimated marginal means), now available on CRAN.
Note that:

  * Really old R scripts that used **lsmeans** will still work with **emmeans** after making 
    minor changes (use `emmeans:::convert_scripts()`). 
  * Existing objects created with really old versions of **lsmeans** can be converted to work 
    with the new package via `emmeans:::convert_workspace()`. 
  * See `vignette("transition-from-lsmeans", "emmeans")` for more details.
  * Datasets `auto.noise`, etc. are now provided only in the **emmeans** package.