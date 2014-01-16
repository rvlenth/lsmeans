# Test cases for new lsmeans design

# subset of warpbreaks
ws = seq(2,53,by=3)

# add a fake covariate
mywarp = transform(warpbreaks, x = rnorm(54))

### stats classes ...
warp.lm = lm(breaks ~ poly(x,3) + wool*tension, data = mywarp, subset = ws)
warp.log = lm(log(breaks) ~ poly(x,3) + wool*tension, data = mywarp, subset=ws)
warp.with = with(mywarp, lm(breaks ~ poly(x,3) + wool*tension, subset = ws))

# A case with empty cells
AM.BL = c(10:18,28:36) # two cells of warpbreaks data
warp.sing = lm(breaks ~ wool*tension, data = warpbreaks, subset = -AM.BL)

data(Oats, package="nlme")
Oats.aov = aov(yield ~ Variety*factor(nitro), data=Oats)

# will fail
Oats.aove = aov(yield ~ Variety + factor(nitro) + Error(Block/Variety), data=Oats)

### nlme classes ...
library(nlme)
Oats.lme = lme(yield ~ Variety*factor(nitro), ~1|Block/Variety, 
               data=Oats)
Oats.lmep = lme(yield ~ Variety*poly(nitro,2), ~1|Block/Variety, 
                data=Oats)
warp.gls = gls(breaks ~ x + I(x^2) + I(x^3) + wool*tension, subset = ws,
                data = mywarp, correlation = corAR1())
# Following will fail because lme objects don't contain coefs info for poly())
warp.gls2 = gls(breaks ~ poly(x,3) + wool*tension, subset = ws,
               data = mywarp, correlation = corAR1())

# lme4
library(lme4)
Oats.lmer = lmer(yield ~ Variety*factor(nitro) + (1|Block/Variety), 
                data=Oats)
Oats.lmerp = lmer(yield ~ Variety*poly(nitro,2) + (1|Block/Variety), 
                data=Oats)
chick.Time = lmer(weight ~ Time*Diet + (0+Time|Chick), data=ChickWeight)
chick.T.origin = lmer(weight ~ Time + Time:Diet + (0+Time|Chick), data=ChickWeight)
chick.logTime = lmer(weight ~ log(Time+1)*Diet + (0+log(Time+1)|Chick), data=ChickWeight)

# MASS
library(MASS)
warp.rlm = rlm(breaks ~ poly(x,3) + wool*tension, data = mywarp)
warp.lqs = lqs(breaks ~ poly(x,3) + wool*tension, data = mywarp)


# Multivariate version of Oats
Oats.mult = with(Oats, expand.grid(Variety=unique(Variety), Block=unique(Block)))
Oats.mult$yield = matrix(Oats$yield, ncol=4, byrow=TRUE)
Oats.mult.lm = lm(yield ~ Block + Variety, data = Oats.mult)
                 

# Load test files
library(plyr)
source("R/recover.data.R")
source("R/ref.grid.R")
source("R/KRstuff.R")
source("R/lsm.basis.R")
source("R/lsm-contr.R")
source("R/lsmeans.R")

# YET TO DO...
# ***DONE*** Warning for effect in an interaction
# *** DECIDED ITS GOOD AS-IS *** Names of list results
# *** DONE *** By variables in contrasts
# *** DONE *** New print method to separate by 'by'?
# *** DONE *** Annotations showing adjust procedure, etc.
# *** DONE *** cld method(s)
# Adjustments to CI methods?
# *** DONE *** Trends -- lstrends function
# Revise NAMESPACE to export only important stuff

