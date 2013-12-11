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
warp.gls = gls(breaks ~ poly(x,3) + wool*tension, subset = ws,
               data = mywarp, correlation = corAR1())

# lme4
library(lme4)
Oats.lmer = lmer(yield ~ Variety*factor(nitro) + (1|Block/Variety), 
                subset = 1:26, data=Oats)
Oats.lmerp = lmer(yield ~ Variety*poly(nitro,2) + (1|Block/Variety), 
                subset = 1:26, data=Oats)