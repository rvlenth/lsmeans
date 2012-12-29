pkgname <- "lsmeans"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('lsmeans')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("fiber")
### * fiber

flush(stderr()); flush(stdout())

### Name: fiber
### Title: Fiber data
### Aliases: fiber
### Keywords: datasets

### ** Examples

require(lsmeans)
fiber.lm <- lm(strength ~ diameter + machine, data=fiber)
lsmeans(fiber.lm, pairwise ~ machine)



cleanEx()
nameEx("lsmeans")
### * lsmeans

flush(stderr()); flush(stdout())

### Name: lsmeans
### Title: Least-squares means
### Aliases: lsmeans print.lsm print.data.frame.lsm lsm glht.lsmlf
### Keywords: models regression htest

### ** Examples

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


### Unbalanced split-plot example ###
#-- The imbalance is imposed deliberately to illustrate that
#-- the variance estimates become biased
require(nlme)
Oats.lme <- lme(yield ~ factor(nitro) + Variety, random = ~1 | Block/Variety, 
    subset = -c(1,2,3,5,8,13,21,34,55), data=Oats)
lsmeans(Oats.lme, list(poly ~ nitro, pairwise ~ Variety))

# Compare with lmer result (if pbkrtest installed, provides df, bias-adjusted SEs)
if (require(pbkrtest)) {
  require(lme4)
  Oats.lmer <- lmer(yield ~ factor(nitro) + Variety + (1 | Block/Variety), 
                   subset = -c(1,2,3,5,8,13,21,34,55), data=Oats)
  lsmeans(Oats.lmer, list(poly ~ nitro, pairwise ~ Variety))
}

# Using in conjunction with 'glht' (note -- this does not use adjusted vcov)
if (require(multcomp)) {
  # calling 'glht' from 'lsmeans' ...
  lsmeans(Oats.lmer, pairwise ~ Variety, glhargs=list(df=9))
  
  # calling 'lsmeans' from 'glht' to get simultaneous CIs
  confint(glht(Oats.lmer, linfct = lsm(~ Variety), df=9))
}

# Custom contrasts
lsmeans(Oats.lmer, my.own ~ Variety, 
  contr = list(my.own = list(G.vs.M = c(1,-1,0), GM.vs.V = c(.5,.5,-1))))




cleanEx()
nameEx("nutrition")
### * nutrition

flush(stderr()); flush(stdout())

### Name: nutrition
### Title: Nutrition data
### Aliases: nutrition
### Keywords: datasets

### ** Examples

require(lsmeans)
nutr.aov <- aov(gain ~ (group + age + race)^2, data = nutrition)

# Summarize predictions for age group 3
nutr.lsm <- lsmeans(nutr.aov, list(pairwise ~ group|race, pairwise ~ race|group),
                   at = list(age="3"))
                   
with(nutr.lsm[[1]], interaction.plot(group, race, lsmean, type="b"))

# Hispanics seem exceptional; but, this doesn't test out due to very sparse data
print(nutr.lsm, omit=3)



cleanEx()
nameEx("pairwise.lsmc")
### * pairwise.lsmc

flush(stderr()); flush(stdout())

### Name: pairwise.lsmc
### Title: Contrast families
### Aliases: pairwise.lsmc revpairwise.lsmc poly.lsmc trt.vs.ctrl.lsmc
###   trt.vs.ctrl1.lsmc trt.vs.ctrlk.lsmc
### Keywords: models regression htest

### ** Examples

### View orthogonal polynomials for 4 levels
poly.lsmc(1:4)

## Not run: 
##D ### Setting up a custom contrast function
##D helmert.lsmc <- function(levs, ...) {
##D   M <- as.data.frame(contr.helmert(levs))
##D   names(M) <- paste(levs[-1],"vs earlier")
##D   attr(M, "desc") <- "Helmert contrasts"
##D   M
##D }
##D lsmeans(Oats.lme, helmert ~ Variety)
## End(Not run)




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
