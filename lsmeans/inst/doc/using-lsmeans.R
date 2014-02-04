### R code from vignette source 'using-lsmeans.rnw'

###################################################
### code chunk number 1: using-lsmeans.rnw:27-31
###################################################
library("lsmeans")
#library("multcomp")
library("lme4")
options(width=88, digits=5, prompt="R> ", cont="R+ ")


###################################################
### code chunk number 2: kbdscat
###################################################
typing <- data.frame(
  keybd = rep(c("A","B","C"), each=4),
  hours = c(60,72,61,50, 54,68,66,59, 56,56,55,51),
  pain =  c(85,95,69,58, 41,74,71,52, 41,34,50,40))
library("lattice")
xyplot(pain ~ hours | keybd, data = typing, layout = c(3,1))


###################################################
### code chunk number 3: using-lsmeans.rnw:141-142
###################################################
typing.lm <- lm(pain ~ hours + keybd, data = typing)


###################################################
### code chunk number 4: using-lsmeans.rnw:145-146
###################################################
( typing.rg <- ref.grid(typing.lm) )


###################################################
### code chunk number 5: using-lsmeans.rnw:149-150
###################################################
summary(typing.rg)


###################################################
### code chunk number 6: using-lsmeans.rnw:153-154
###################################################
( typing.lsm <- lsmeans(typing.rg, "keybd") )


###################################################
### code chunk number 7: using-lsmeans.rnw:159-160
###################################################
ref.grid(typing.lm, cov.reduce = median)


###################################################
### code chunk number 8: using-lsmeans.rnw:163-165
###################################################
typing.rg2 <- ref.grid(typing.lm, at = list(hours = c(50,60)))
lsmeans(typing.rg2, c("keybd","hours"))


###################################################
### code chunk number 9: using-lsmeans.rnw:168-170
###################################################
lsmeans(typing.rg2, "keybd")
lsmeans(typing.rg2, "hours")


###################################################
### code chunk number 10: using-lsmeans.rnw:176-177
###################################################
( typing.pairs <- pairs(typing.lsm) )


###################################################
### code chunk number 11: using-lsmeans.rnw:180-181
###################################################
cld(typing.lsm, alpha = .10)


###################################################
### code chunk number 12: using-lsmeans.rnw:186-187
###################################################
contrast(typing.lsm, "eff")


###################################################
### code chunk number 13: using-lsmeans.rnw:193-194
###################################################
summary(typing.pairs, adjust = "none")


###################################################
### code chunk number 14: using-lsmeans.rnw:203-206
###################################################
library("multcomp")
typing.glht <- glht(typing.lm, typing.pairs)
summary(typing.glht)


###################################################
### code chunk number 15: typing-glht-plot
###################################################
plot(typing.glht)


###################################################
### code chunk number 16: using-lsmeans.rnw:218-219
###################################################
confint(glht(typing.lm, lsm("keybd")))


###################################################
### code chunk number 17: using-lsmeans.rnw:222-223 (eval = FALSE)
###################################################
## lsm("keybd", contr="pairwise")


###################################################
### code chunk number 18: using-lsmeans.rnw:230-232
###################################################
lsmeans(typing.lm, specs = "keybd", by = "hours", 
    at = list(hours = c(50, 60)), contr = "trt.vs.ctrl1")


###################################################
### code chunk number 19: using-lsmeans.rnw:237-245 (eval = FALSE)
###################################################
## lsmeans(typing.lm, specs = ~ keybd, by = "hours", 
##     at = list(hours = c(50, 60)), contr = "trt.vs.ctrl1")
##     
## lsmeans(typing.lm, specs =  ~ keybd | hours, 
##     at = list(hours = c(50, 60)), contr = "trt.vs.ctrl1")
##     
## lsmeans(typing.lm, specs =  trt.vs.ctrl1 ~ keybd | hours, 
##     at = list(hours = c(50, 60)))


###################################################
### code chunk number 20: using-lsmeans.rnw:252-254
###################################################
noise.lm <- lm(noise ~ size*type*side, data = auto.noise)
anova(noise.lm)


###################################################
### code chunk number 21: using-lsmeans.rnw:257-258
###################################################
ref.grid(noise.lm)


###################################################
### code chunk number 22: noise-ip
###################################################
lsmip(noise.lm, size ~ type | side)


###################################################
### code chunk number 23: using-lsmeans.rnw:269-272 (eval = FALSE)
###################################################
## lsmip(noise.lm, size ~ type * side)  # 1 panel, 3 curves, 2*2 = 4 x values
## lsmip(noise.lm, type * side ~ size)  # 1 panel, 2*2 = 4 curves, 3 x values
## lsmip(noise.lm, type ~ side | size)  # 3 panels, 2 curves, 2 x values


###################################################
### code chunk number 24: using-lsmeans.rnw:276-277
###################################################
lsmeans(noise.lm, pairwise ~ type)


###################################################
### code chunk number 25: using-lsmeans.rnw:286-287
###################################################
with(auto.noise, tapply(noise, type, mean))


###################################################
### code chunk number 26: using-lsmeans.rnw:292-293
###################################################
lsmeans(noise.lm, pairwise ~ type | size*side)[[2]]


###################################################
### code chunk number 27: oats-intplot
###################################################
data(Oats, package = "nlme")
library("lme4")
Oats.lmer <- lmer(yield ~ Variety*factor(nitro) + (1|Block/Variety), 
    data = Oats)
anova(Oats.lmer)
lsmip(Oats.lmer, Variety ~ nitro)


###################################################
### code chunk number 28: using-lsmeans.rnw:316-319
###################################################
Oats.add <- lmer(yield ~ Variety + factor(nitro) + (1|Block/Variety),
    data = Oats)
lsmeans(Oats.add, list(revpairwise ~ Variety,  poly ~ nitro))


###################################################
### code chunk number 29: using-lsmeans.rnw:322-324
###################################################
Oats.poly <- lmer(yield ~ Variety + poly(nitro, 2) + (1 | Block/Variety), 
    data=Oats)


###################################################
### code chunk number 30: using-lsmeans.rnw:327-330
###################################################
Oats.poly.rg <- ref.grid(Oats.poly, at = list(nitro = c(0, .2, .4, .6)))
lsmeans(Oats.poly.rg, ~ Variety)
lsmeans(Oats.poly.rg, ~ nitro)


###################################################
### code chunk number 31: oatspoly-intplot
###################################################
lsmip(Oats.poly, Variety ~ nitro, cov.reduce = FALSE)


###################################################
### code chunk number 32: nutr-intplot
###################################################
nutr.lm <- lm(gain ~ (age + group + race)^2, data = nutrition)
lsmip(nutr.lm, race ~ age | group)
lsmeans(nutr.lm, ~ group*race)


###################################################
### code chunk number 33: using-lsmeans.rnw:365-366
###################################################
lsmeans(nutr.lm, ~ group*race, at = list(age = "3"))


###################################################
### code chunk number 34: using-lsmeans.rnw:370-373
###################################################
wtavg <- function(coefs, lev) (23*coefs[1,] + 64*coefs[2,])/87
nutr.lsm <- lsmeans(nutr.lm, ~ group * race, fac.reduce = wtavg,
    at = list(age=c("2","3"), race=c("Black","White")))


###################################################
### code chunk number 35: using-lsmeans.rnw:376-379
###################################################
nutr.lsm    
pairs(nutr.lsm, by = "race")
pairs(nutr.lsm, by = "group")


###################################################
### code chunk number 36: chick-plot
###################################################
xyplot(weight~Time | Diet, groups = ~ Chick, data=ChickWeight, type="o", 
       layout=c(4,1))


###################################################
### code chunk number 37: using-lsmeans.rnw:399-401
###################################################
Chick.lmer <- lmer(weight ~ Diet * Time + (0 + Time | Chick), 
    data = ChickWeight)


###################################################
### code chunk number 38: using-lsmeans.rnw:404-405
###################################################
cld (lstrends (Chick.lmer, ~ Diet, var = "Time"))


###################################################
### code chunk number 39: using-lsmeans.rnw:410-413
###################################################
Chick.lmer2 <- lmer(weight ~ Diet * log(Time + 1) + 
    (0 + log(Time + 1) | Chick),  data = ChickWeight)
cld (lstrends (Chick.lmer2, ~ Diet, var = "log(Time + 1)"))


###################################################
### code chunk number 40: using-lsmeans.rnw:416-417
###################################################
cld (lstrends (Chick.lmer2, ~ Diet, var = "Time"))


###################################################
### code chunk number 41: using-lsmeans.rnw:424-425
###################################################
MOats.mlm <- lm(yield ~ Block + Variety, data = MOats)


###################################################
### code chunk number 42: using-lsmeans.rnw:429-430
###################################################
ref.grid(MOats.mlm)


###################################################
### code chunk number 43: using-lsmeans.rnw:433-435
###################################################
MOats.rg <- ref.grid(MOats.mlm, mult.levs = list(nitro = c(0,.2,.4,.6)))
MOats.rg


###################################################
### code chunk number 44: using-lsmeans.rnw:438-440
###################################################
( MOats.lsm <- lsmeans(MOats.rg, ~ nitro | Variety) )
( MOats.pcon <- contrast(MOats.lsm, "poly") )


###################################################
### code chunk number 45: using-lsmeans.rnw:443-444
###################################################
pairs(MOats.pcon, by = "contrast")


###################################################
### code chunk number 46: using-lsmeans.rnw:451-454
###################################################
cbpp$obs <- 1:nrow(cbpp)
cbpp.glmer <- glmer(cbind(incidence, size - incidence) 
   ~ period + (1 | herd) + (1 | obs),  family = binomial,  data = cbpp)


###################################################
### code chunk number 47: using-lsmeans.rnw:457-464
###################################################
cbpp.lsm <- lsmeans(cbpp.glmer, ~ period)
cbpp.sum <- summary(cbpp.lsm)
cbpp.sum$pred.incidence <- 1 - 1 / (1 + exp(cbpp.sum$lsmean))
cbpp.sum
cbpp.con <- summary(contrast(cbpp.lsm, "trt.vs.ctrl1"))
cbpp.con$odds.ratio <- exp(cbpp.con$estimate)
cbpp.con


