R package **lsmeans**: Least-squares means (predicted marginal means)
====

[![Build Status](https://travis-ci.org/rvlenth/lsmeans.svg?branch=master)](https://travis-ci.org/rvlenth/lsmeans)
[![cran version](http://www.r-pkg.org/badges/version/lsmeans)](https://cran.r-project.org/package=lsmeans)
[![downloads](http://cranlogs.r-pkg.org/badges/lsmeans)](http://cranlogs.r-pkg.org/badges/lsmeans)
[![total downloads](http://cranlogs.r-pkg.org/badges/grand-total/lsmeans)](http://cranlogs.r-pkg.org/badges/grand-total/lsmeans)
[![Research software impact](http://depsy.org/api/package/cran/lsmeans/badge.svg)](http://depsy.org/package/r/lsmeans)

## Features
* For an overview, see the [vignette](https://CRAN.R-project.org/package=lsmeans/vignettes/using-lsmeans.pdf) -- also available via `vignette("using-lsmeans", package = "lsmeans")`
* Least-squares means (a.k.a. predicted marginal means) are derived by using a model to make predictions over a regular grid of pridictor combinations (called a *reference grid*). These predictions may possibly be averaged (typically with equal weights) over one or more of the predictors. Such marginally-averaged predictions are useful for describing the results of fitting a model, particularly in presenting the effects of factors. The **lsmeans** package can easily produce these results, as well as various graphs of them (interaction-style plots and side-by-side intervals).
* Estimation and testing of pairwise comparisons of LS means, and several other types of contrasts, are provided. There is also a `cld` method for display of grouping symbols.
* Two-way support of the `glht` function in the **multcomp** package.
* For models where continuous predictors interact with factors, the package's `lstrends` function works in terms of a reference grid of predicted slopes of trend lines for each factor combination.
* Incorporates support for many types of models, including those in **stats** package (`lm`, `glm`, `aov`, `aovlist`), linear and genearalized linear mixed models (e.g., **nlme**, **lme4**, **afex**), ordinal-response models (e.g., **ordinal**, **MASS**), survival analysis (e.g., **survival**, **coxme**), generalized estimating equations (**gee**, **geepack**), and others. See `help("models"", package="lsmeans")`
* Various Bayesian models (**carBayes**, **MCMCglmm**, **MCMCpack**) are supported by way of creating a posterior sample of least-squares means or contrasts thereof, which may then be examined using tools such as in the **coda** package.
* Package developers may provide **lsmeans** support for their models by providing `recover.data` and `lsm.basis` methods. See `vignette("extending", package = "lsmeans")`

## Installation
* To install latest version from CRAN, run 
```
install.packages("lsmeans")
```
Release notes for the latest CRAN version are found at [https://CRAN.R-project.org/package=lsmeans/NEWS](https://CRAN.R-project.org/package=lsmeans/NEWS) -- or do `news(package = "lsmeans")` for notes on the version you have installed.

* To install the latest development version from Github, have the newest **devtools** package installed, then run
```
devtools::install_github("rvlenth/lsmeans", dependencies = TRUE)
```
For latest release notes on this development version, see the [NEWS file](https://github.com/rvlenth/lsmeans/blob/master/inst/NEWS)
