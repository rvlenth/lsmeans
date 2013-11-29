\name{lsmeans}
\alias{lsmeans}
\alias{print.lsm}
\alias{print.data.frame.lsm}
\alias{lsm}
\alias{glht.lsmlf}

\title{Least-squares means}
\description{
Compute least-squares means for specified factors or factor combinations in a linear model,
and optionally comparisons or contrasts among them.
}
\usage{
lsmeans(object, specs, adjust = c("auto","tukey","sidak","scheffe", p.adjust.methods),
  conf = .95, at, trend, contr = list(), 
  cov.reduce = function(x, name) mean(x), 
  fac.reduce = function(coefs, lev) apply(coefs, 2, mean), 
  glhargs = NULL, lf = FALSE, mlf = rep(1,nresp)/nresp, ...)
  
\method{print}{lsm}(x, omit = NULL, ...)
\method{print}{data.frame.lsm}(x, ...)

lsm(...)
\method{glht}{lsmlf}(model, linfct, ...)
}
\arguments{
  \item{object}{
A \code{lm}, \code{aov} (with no \code{Error} component), \code{glm}, \code{lme}, \code{gls}, \code{lmer}, or \code{glmer} object having at least one fixed factor among the predictors.
}
  \item{specs}{
A formula, or a list of formulas, specifying the desired families of least-squares means.
The right-hand side of each formula specifies the desired predictors for which we want least-squares means. The optional left-hand side specifies
what kind of comparisons or contrasts are desired. For example, \code{~ treatment} requests least-squares means for each level of \code{treatment}, and \code{pairwise ~ treatments} requests those results, plus pairwise comparisons among them.
As another example, in a three-factor model, \code{trt.vs.ctrl1 ~ A | B:C} requests least-squares means for all combinations of factors \code{A}, \code{B}, and \code{C}, as well as treatment-minus-control comparisons of \code{A} for each combination of \code{B} and \code{C}, where the first level of \code{A} is considered the control level.
}
  \item{adjust}{
Adjustment method for the p values of tests of contrasts.
\code{"auto"} uses the method returned in the \code{"adjust"} attribute of the contrast function;
\code{"tukey"} computes p values using the Studentized range distribution with the number of means in the family;
\code{"sidak"} replaces each p value by \code{1 - (1 - p)^c}, where c is the number of contrasts;
\code{"scheffe"} computes p values from the F distribution, according to the Scheffe criterion;
\code{"none"} makes no adjustments to the p values.
In addition, the Bonferroni-inequality-based adjustment methods in \code{\link{p.adjust}} are available; many of those involve step-down testing with a stopping rule, whereas \code{"tukey"}, \code{"sidak"}, \code{"bonferroni"}, and \code{"none"} anticipate that all contrasts will be tested.
In many cases, these adjustments are only approximate, especially when the degrees of freedom vary
greatly within the family of comparisons. For more accurate adjustments, use \code{glhargs} instead.
}
  \item{conf}{
Desired confidence level for intervals. For robustness, you may specify either a fraction or a percentage; i.e., \code{.975} and \code{97.5} yield the same results.
}
  \item{at}{
An optional named list or named vector of predictor values for the reference grid. If no value is found in \code{at} for a particular predictor, then if it is a factor, all its levels are used, and if a covariate, we use the results of \code{cov.reduce}.  
}
  \item{trend}{
If non-missing and \code{trend} is the name of a numeric predictor or a term in \code{object}'s model,
then the linear trends for this predictor are estimated and compared instead of the least-squares means. This is useful in models where the specified variable or term interacts with one or more of the factors in the right-hand side of \code{specs}. In multi-factor cases, sets of 
regression coefficients are averaged using \code{fac.reduce}, just as with least-squares means.
}
  \item{contr}{
An optional named \code{list}. Each entry is itself a \code{list} or a \code{data.frame} specifying contrast coefficients. If the left-hand side of a formula in \code{specs} matches a name in \code{contr}, then those contrasts are estimated with the specified least-squares means. An error will result if the length or one or more contrast vectors mismatches the number of levels of the factor or factor combination. Actually, it is not necessary that they be contrasts; one may use this argument to estimate arbitrary linear combinations of the least-squares means.
}
  \item{cov.reduce}{
A logical value, or a function with arguments \code{x} and \code{name} that should return the value to use in prediction for the covariate with name \code{name} and values \code{x}. By default, the mean is used. The \code{name} argument may help distinguish one covariate from another. 

If \code{cov.reduce} is logical, a value of \code{TRUE} is the same as the default, and a value of \code{FALSE} calls \code{\link{unique}} for each covariate (be careful using this if covariates have more than a few levels).
}
  \item{fac.reduce}{
A function of \code{coefs} and \code{lev} where \code{lev} is the level of a factor or factor combination at which a least-squares mean is calculated.
The argument \code{coefs} is a matrix whose rows correspond to the combinations of all factors in the model other than those involved in the \code{lsmeans} specification. Each row has the coefficients for the linear combination of the regression coefficients to be used in that case. By default, these rows are averaged together (mimicking SAS), but the user may override that behavior. Besides \code{lev}, \code{names(lev)} will provide the name of the factor or factor combination, and the \code{row.names} of \code{coefs} provide the levels of the extraneous factors.
}
  \item{glhargs}{
If this is a \code{list}, the object and specified contrasts are passed to the function \code{\link[multcomp]{glht}} in the \pkg{multcomp} package, with the contents of \code{glhargs} as additional arguments. (If you do not wish to provide additional arguments, use \code{glsargs=list()}.) If \code{glhargs} is left at \code{NULL}, 
%%%or if the \pkg{multcomp} package is not installed, 
then \code{glht} is not called, and the contrast results are produced internally by \code{lsmeans}. This argument affects only the results from contrasts and not those for the lsmeans themselves. Note: If \code{glhargs} is used, the \code{adjust} argument is ignored. 
If only one set of lsmeans or contrasts thereof is desired, an alternative is to call \code{glht} with an \code{lsm} call as the \code{linfct} argument.
}
\item{lf}{
If \code{TRUE}, the linear functions of the regression coefficients are returned unevaluated.
}
\item{model}{
\code{model} argument for \code{\link[multcomp]{glht}}
}
\item{linfct}{
\code{linfct} argument for \code{\link[multcomp]{glht}}
}
\item{mlf}{
   Numeric vector defining a linear function of a multivariate response (i.e., \code{y \%*\% mlf}
   where \code{y} is the response variable of dimension \code{nresp}). The default is the mean.)
}
  \item{\dots}{
Additional argument(s) passed to the contrast function, or to \code{\link[multcomp]{glht}}; see Details.
}
\item{x}{Object of class \code{"lsm"}}
\item{omit}{Indexes of elements of \code{x} that you do not want printed.}
}
\details{
Least-squares means, popularized by SAS, are predictions from a linear model at combinations
of specified factors. SAS's documentation describes them as ``predicted population margins---that is, they estimate the marginal means over a balanced population'' (SAS Institute 2012). In generalized linear models, least-squares means are marginal linear predictions that can be transformed back to the response scale via the inverse-link function.
Unspecified factors and covariates are handled by summarizing the predictions
over those factors and variables. For example, if the fitted model has formula \code{response ~ x1 + x2 + treat}
where \code{a1} and \code{x2} are numeric and \code{treat} is a factor, the least-squares means will be the predicted response for each treatment, at some specified values of \code{x1} and \code{x2}. By default, the means of the two covariates will be used, resulting in what ANOVA textbooks often call the adjusted means. We may use that \code{at} argument to instead make predictions at other values of \code{x1} and \code{x2}.

Now consider the model \code{response ~ A + B + A:B}, where \code{A} and \code{B} are both factors. If we ask for least-squares means for \code{A}, then at each level of \code{A} we are faced with a different prediction for each level of \code{B}. Blind (and default) use of least-squares means would result in these predictions being averaged together with equal weight, and this may be inappropriate, especially when the interaction effect is strong. (A warning is generated in such potentially inappropriate cases.) Like most statistical calculations, it is possible to use least-squares means inappropriately. The \code{fac.reduce} argument at least expands one's options in producing meaningful results in multi-factor situations. 

If you give multiple values for a covariate in \code{at} or as a result of \code{cov.reduce}, the covariate is treated just like a factor, in that predictions are made at each combination of its values with other factors and multi-valued covariates, and averaged as appropriate using \code{fac.reduce}. If you want to see separate predictions for each value, include the covariate in \code{specs}. (This differs from what SAS does with non-unique \code{at} values.) 

One other note concerning covariates: One must be careful with covariates that depend on one another. For example, if a model contains covariates \code{x} and \code{xsq} where \code{xsq = x^2}, the default behavior will make predictions at \code{x = mean(x)} and \code{xsq = mean(xsq)}, which probably isn't a valid combination (we need \code{x = mean(x)} and \code{xsq = mean(x)^2}). The inconsistency is avoided if the model specifis \code{poly(x,2)} (or even \code{x + I(x^2)}) instead of \code{x + xsq}, because then only \code{x} is averaged, and everything remains consistent.

The built-in contrast methods that can be used in \code{specs} formulas are \code{pairwise}, \code{revpairwise}, \code{poly}, \code{trt.vs.ctrl}, \code{trt.vs.ctrl1}, and \code{trt.vs.ctrlk}. They are implemented as functions having the same names with \code{.lsmc} added (e.g., \code{\link{pairwise.lsmc}}). Users may write additional \code{.lsmc} functions that generate custom families of contrasts. See the documentation for \code{\link{pairwise.lsmc}} for an example.

Experimental feature: If a contrast family named \code{cld} is specified, then it is the same as \code{pairwise}, with the side effect that a column named \code{.group} is added to the table of least-squares means. This column is a compact letter display (CLD): means associated with one or more of the same letters have P values greater than \code{1 - conf} in the selected pairwise-comparison tests. For more details on CLDs, see the documentation for \code{\link[multcomp]{cld}} in the \pkg{multcomp} package.

When \code{trend} is specified and it names a predictor variable in the model, the trend is estimated using a difference quotient over 1/1000 the range of the predictor (centered at the predictor's \code{at} value or its \code{cov.reduce} result). If \code{trend} instead names a model term, then the trend is determined symbolically. It could be informative to run with \code{lf=TRUE}. See the examples.

Degrees of freedom are currently not provided for \code{lme} or \code{glme} objects, or for \code{mer} objects arising from generalized linear models; in those cases, asymptotic results are printed, and this fact is emphasized by displaying \code{NA} for the defrees of freedom. For linear \code{mer} objects, degrees of freedom are computed using the Kenward and Roger (1997) method, provided the \code{pbkrtest} package is installed (the package is loaded if needed). Moreover, in that case, the adjusted covariance matrix from the \code{vcovAdj()} function in the \code{pbkrtest} package is used to calculate standard errors. See Halekoh and H?jsgaard (2012) and the documentation for \code{\link[pbkrtest]{KRmodcomp}} for more details. Degrees of freedom are not known to \code{\link[multcomp]{glht}} except in the case of \code{lm} objects. Note that \code{glht} requires the degrees of freedom to be an integer; accordingly, if \code{df} is included in \code{glhargs}, it is coerced to an integer value of \code{max(1, as.integer(df + .2))}. That is, it is forced to be at least 1, and it is taken up to the next higher integer if its fractional part is at least .8.

For users of the package \pkg{multcomp}, linear functions corresponding to \code{specs} may be passed to \code{\link[multcomp]{glht}} in two different ways. One is to call \code{glht} from within \code{lsmeans} by specifying \code{glhtargs}. The other way is to call \code{lsmeans} from within \code{glht} by supplying a call to \code{lsm} in its \code{linfct} argument. The first way provides \code{summary.glht} objects for all contrast results (but not for the lsmeans). The second way passes \emph{just one} set of linear functions to \code{glht}; in particular, only the first element of \code{specs} is used, and the last set of linear functions from that specification are given to \code{glht}. Thus, \code{glht(model, lsm(~treat))} will estimate the lsmeans of treatments, while \code{glht(model, lsm(pairwise~treat))} will estimate the pairwise comparisons.

If the model contains a matrix among its predictors, each column is averaged using the function specified in \code{cov.reduce}. There is no provision for matrices in the \code{at} argument.
}
\value{
If \code{lf==FALSE}, \code{lsmeans} returns an object of class \code{"lsm"}, which inherits from \code{"list"}. Each element of the list is either of class \code{data.frame.lsm} or \code{"summary.glht"} (see the documentation for \code{\link[multcomp]{glht}}). (The latter occur only if \code{glhargs} is non-NULL.) Each element summarizes a family of least-squares means or contrasts among them. Each \code{data.frame.lsm} element is an extension of a \code{\link{data.frame}} and contains lsmeans or contrast estimates and associated quantities.

If \code{lf==TRUE}, \code{lsmeans} returns a list of matrices containing the linear functions generated by \code{specs}. 

The \code{"lsm"} and \code{"data.frame.lsm"} classes each have only a \code{print} method, which exist to display results in the desired format.

The \code{lsm} function returns an object of class \code{"lsmlf"}, and may be used in a manner similar to \code{\link[multcomp]{mcp}} in a call to \code{\link[multcomp]{glht}}. This is implemented via the provided S3 method for \code{glht} for \code{"lsmlf"} objects. See the examples.
}

%%%%\note{}

\references{
Halekoh, U. and H?jsgaard, S. (2012)
A Kenward-Roger Approximation and parametric bootsrap methods for tests in linear mixed models -- the R package \code{pbkrtest}, submitted. %%%%%\emph{Journal of Statistical Software}.

Kenward, M.G. and Roger, J.H. (1997)
Small sample inference for fixed effects from restricted maximum likelihood,
\emph{Biometrics}, 53, 983--997.

SAS Institute Inc. (2012)
Online documentation; Shared concepts; LSMEANS statement,
\url{http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_introcom_a0000003362.htm}, accessed August 15, 2012.
}

\author{
Russell V. Lenth, The University of Iowa
}

\seealso{
\code{\link{pairwise.lsmc}}, \code{\link[multcomp]{glht}}
}

\examples{
require(lsmeans)

### Covariance example (from Montgomery Design (8th ed.), p.656)
# Uses supplied dataset 'fiber'
fiber.lm <- lm(strength ~ diameter + machine, data = fiber)

# adjusted means and comparisons, treating machine C as control
lsmeans (fiber.lm, trt.vs.ctrlk ~ machine)


### Factorial experiment
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#-- We only need to see the wool*tension means listed once ...
print(lsmeans (warp.lm,  list(pairwise ~ wool | tension,  poly ~ tension | wool)),
    omit=3)

# Compact letter display
lsmeans(warp.lm, cld ~ tension | wool)


### Unbalanced split-plot example ###
#-- The imbalance is imposed deliberately to illustrate that
#-- the variance estimates become biased
require(nlme)
Oats.lme <- lme(yield ~ factor(nitro) + Variety, random = ~1 | Block/Variety, 
    subset = -c(1,2,3,5,8,13,21,34,55), data=Oats)
lsmeans(Oats.lme, list(poly ~ nitro, pairwise ~ Variety))

# Model with a quadratic trend for 'nitro'
Oatsq.lme <- update(Oats.lme, . ~ nitro + I(nitro^2) + Variety)

# Predictions at each unique 'nitro' value in the dataset
lsmeans(Oatsq.lme, ~ nitro, cov.reduce = FALSE)

# lsmeans for 'Variety' at the average 'nitro' value.
lsmeans(Oatsq.lme, ~ Variety)
# These differ quite a bit from those obtained from 'Oats.lme', because each
# is a single prediction rather than the average of 4 predictions

# lsmeans for 'Variety', using the four unique 'nitro' values
lsmeans(Oatsq.lme, ~ Variety, cov.reduce = FALSE)
# These are pretty close to those from 'Oats.lme' because both are based
# on averaging predictions for the 12 


# Compare with lmer result 
if (require(lme4)) {
    
    Oats.lmer <- lmer(yield ~ factor(nitro) + Variety + (1 | Block/Variety), 
                      subset = -c(1,2,3,5,8,13,21,34,55), data=Oats)
    # If pbkrtest installed, provides df, bias-adjusted SEs
    if (require(pbkrtest)) {
        lsmeans(Oats.lmer, list(poly ~ nitro, pairwise ~ Variety))
    }
    
    # Using in conjunction with 'glht' (note -- this does not use adjusted vcov)
    # calling 'glht' from 'lsmeans' ...
    lsmeans(Oats.lmer, pairwise ~ Variety, glhargs=list(df=9))

    # calling 'lsmeans' from 'glht' to get simultaneous CIs
    confint(glht(Oats.lmer, linfct = lsm(~ Variety), df=9))

    # Custom contrasts
    lsmeans(Oats.lmer, my.own ~ Variety, 
            contr = list(my.own = list(G.vs.M = c(1,-1,0), GM.vs.V = c(.5,.5,-1))))
}
    

#### Examples with trends...
# Model with interaction
fiber.lm2 = lm(strength ~ diameter * machine, data = fiber)

# Compare the linear trends for diameter
lsmeans(fiber.lm2, pairwise ~ machine, trend = "diameter")

# Model with log(diameter) as the covariate
fiber.lm3 = lm(strength ~ log(diameter) * machine, data = fiber)
# Compare the fitted linear trends for log(diameter)
lsmeans(fiber.lm3, pairwise ~ machine, trend = "log(diameter)")

# Compare the fitted linear trends for diameter itself 
# - this is done via a diff quotient - compare with fiber.lm2 results
lsmeans(fiber.lm3, pairwise ~ machine, trend = "diameter")

\dontrun{
# Examine the linear functions generated for these examples
lsmeans(fiber.lm2, ~ machine, trend = "diameter", lf = TRUE)
lsmeans(fiber.lm3, ~ machine, trend = "log(diameter)", lf = TRUE)
lsmeans(fiber.lm3, ~ machine, trend = "diameter", lf = TRUE)
}
}
\keyword{ models }
\keyword{ regression }
\keyword{ htest }
