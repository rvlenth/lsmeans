# Test cases for new lsmeans design

# subset of warpbreaks
ws = seq(2,53,by=3)

# add a fake covariate
mywarp = transform(warpbreaks, x = rnorm(54))

### stats classes ...
warp.data = lm(breaks ~ poly(x,3) + wool*tension, data = mywarp, subset = ws)
warp.log = lm(log(breaks) ~ poly(x,3) + wool*tension, data = mywarp, subset=ws)
warp.with = with(mywarp, lm(breaks ~ poly(x,3) + wool*tension, subset = ws))
attach(mywarp)
warp.att = lm(breaks ~ poly(x,3) + wool*tension, subset = ws)
detach()
# A case with empty cells
AM.BL = c(10:18,28:36) # two cells of warpbreaks data
warp.sing = lm(breaks ~ wool*tension, data = warpbreaks, subset = -AM.BL)

data(Oats, package="nlme")
Oats.aov = aov(yield ~ Variety*factor(nitro), subset = 1:26, data=Oats)

# will fail
Oats.aove = aov(yield ~ Variety + factor(nitro) + Error(Block/Variety), subset = 1:24, data=Oats)

### nlme classes ...
library(nlme)
Oats.lme = lme(yield ~ Variety*factor(nitro), ~1|Block/Variety, 
               subset = 1:26, data=Oats)
Oats.lmep = lme(yield ~ Variety*poly(nitro,2), ~1|Block/Variety, 
                subset = 1:26, data=Oats)
warp.gls = gls(breaks ~ x + I(x^2) + I(x^3) + wool*tension, subset = ws,
                data = mywarp, correlation = corAR1())
# Following will fail because lme objects don't contain coefs info for poly())
warp.gls2 = gls(breaks ~ poly(x,3) + wool*tension, subset = ws,
               data = mywarp, correlation = corAR1())

# lme4
library(lme4)
Oats.lmer = lmer(yield ~ Variety*factor(nitro) + (1|Block/Variety), 
                subset = 1:26, data=Oats)
Oats.lmerp = lmer(yield ~ Variety*poly(nitro,2) + (1|Block/Variety), 
                subset = 1:26, data=Oats)

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
# Warning for effect in an interaction
# By variables in contrasts
# New print method to separate by 'by'?
# Annotations showing adjust proceduere, etc.
# cld method(s)
# Adjustments to CI methods?
# Trends
