# Additional tests of lsmeans for various models
# tests1.R tests the features of lsmeans, so we don't
# repeat those here. The primary purpose of these tests
# is to ensure that various models other than lm and its
# relatives are adequately supported. So the main thrust
# to simply confirm that we can build a reference grid
# and extract the needed statistics.

require(lsmeans)

## --- lme and gls ---
require(nlme)

Oats.lme <- lme(yield ~ Variety + factor(nitro), ~1|Block/Variety, 
               data=Oats)
summary(ref.grid(Oats.lme))

warp.gls <- gls(breaks ~ wool*tension, 
                data = warpbreaks, correlation = corAR1())
summary(ref.grid(warp.gls))


## --- lmer, glmer, glmer.nb ---
require(lme4)

Oats.lmer <- lmer(yield ~ Variety + factor(nitro) + (1|Block/Variety), 
               data=Oats)
summary(ref.grid(Oats.lmer))

cbpp.glmer <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             family = binomial, data = cbpp)
summary(ref.grid(cbpp.glmer))

grouse.gnb  <- glmer.nb(TICKS ~ YEAR+HEIGHT+(1|BROOD), 
	data=grouseticks)
summary(ref.grid(grouse.gnb, data = grouseticks))

## --- coxph and survreg ---
require(survival)

#### borrowed from man pages
bladder1 <- bladder[bladder$enum < 5, ] 
blad.cph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
      cluster(id), data = bladder1)
summary(ref.grid(blad.cph))      

lung.sr <- survreg(Surv(time, status) ~ ph.ecog + age + strata(sex), 
	data = lung)      
summary(ref.grid(lung.sr), type = "response")



## --- coxme ---
require(coxme)

eortc.cme <- coxme(Surv(y, uncens) ~ trt + (trt| center) + strata(center), 
	data = eortc)
summary(ref.grid(eortc.cme, at=list(center=c(1:5))))



## --- selected MASS models ---
## Except for polr, these work just because of inheriting from
## other objects
require(MASS)

warp.rlm <- rlm(breaks ~ wool*tension, data = warpbreaks)
summary(ref.grid(warp.rlm))

house.plr <- polr(Sat ~ (Infl + Type + Cont)^2, weights = Freq, 
	data = housing)
summary(ref.grid(house.plr))

quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)
summary(ref.grid(quine.nb1))

bact.pql <- glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
	family = binomial, data = bacteria)
summary(ref.grid(bact.pql), type = "response")
