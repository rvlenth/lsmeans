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
warp.rlm = rlm(breaks ~ poly(x,3) + wool*tension, data = mywarp, subset=ws)
warp.lqs = lqs(breaks ~ poly(x,3) + wool*tension, data = mywarp, subset=ws)


# Multivariate version of Oats
Oats.mult = with(Oats, expand.grid(Variety=unique(Variety), Block=unique(Block)))
Oats.mult$yield = matrix(Oats$yield, ncol=4, byrow=TRUE)
Oats.mult.lm = lm(yield ~ Block + Variety, data = Oats.mult)


# Survival models
library(survival)
Oats.sr = survreg(Surv(sapply(yield,min,100), yield<100) ~ Block + factor(nitro) + Variety, 
    dist="gaussian", data = Oats, subset=13:72)

# CGD data discussed in Therneau & Crowson papaer w/ survival package
cgd.ph <- coxph(Surv(tstart, tstop, status) ~ treat * inherit + 
    sex + age + cluster(id), data = cgd)

cgd.rg <- ref.grid(cgd.ph)
# Here are the "right answers":
predict(cgd.ph, newdata=cgd.rg@grid, se.fit=TRUE)
# Compare with my answers:
summary(cgd.rg)

swiss.lm = lm(Fertility ~ (Agriculture + Education + Examination + 
                    Infant.Mortality)^2, data = swiss)

swiss.rg = ref.grid(swiss.lm, at=list(
    Agriculture = c(25,50,75),
    Examination = c(8,16,24),
    Infant.Mortality = c(15,20,25)))

lsmip(swiss.rg, Agriculture ~ Examination | Infant.Mortality,
      ylab = "Predicted Fertility")

lsmeans(swiss.rg, "Examination")

lsmeans(swiss.rg, "Examination", by = "Infant.Mortality")
