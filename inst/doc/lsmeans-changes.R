### R code from vignette source 'lsmeans-changes.rnw'

###################################################
### code chunk number 1: lsmeans-changes.rnw:35-41
###################################################
library(lsmeans)

### OLD
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.oldlsm <- .old.lsmeans(warp.lm, ~ tension | wool)
class(warp.oldlsm)


###################################################
### code chunk number 2: lsmeans-changes.rnw:43-44
###################################################
class(warp.oldlsm[[1]])


###################################################
### code chunk number 3: lsmeans-changes.rnw:47-50
###################################################
### NEW
warp.lsmobj <- lsmeans(warp.lm, ~ tension | wool)
class(warp.lsmobj)


###################################################
### code chunk number 4: lsmeans-changes.rnw:52-53
###################################################
"lsmobj"


###################################################
### code chunk number 5: lsmeans-changes.rnw:56-57
###################################################
warp.lsmobj


###################################################
### code chunk number 6: lsmeans-changes.rnw:62-64
###################################################
### OLD
warp.oldlsm[[1]]$lsmean


###################################################
### code chunk number 7: lsmeans-changes.rnw:66-69
###################################################
Try <- function(expr) tryCatch(expr, error = function(e) cat("Oops!\n"))
### NEW
Try(warp.lsmobj$lsmean)


###################################################
### code chunk number 8: lsmeans-changes.rnw:72-74
###################################################
### NEW
summary(warp.lsmobj)$lsmean


###################################################
### code chunk number 9: lsmeans-changes.rnw:77-78
###################################################
as.data.frame(summary(warp.lsmobj))


###################################################
### code chunk number 10: lsmeans-changes.rnw:82-86
###################################################
### NEW
warp.l2 <- lsmeans(warp.lm, pairwise ~ tension)
class(warp.l2)
sapply(warp.l2, class)


###################################################
### code chunk number 11: lsmeans-changes.rnw:100-101
###################################################
confint(warp.lsmobj, level = .90)


###################################################
### code chunk number 12: lsmeans-changes.rnw:104-105
###################################################
summary(glht(warp.lm, warp.lsmobj))


###################################################
### code chunk number 13: lsmeans-changes.rnw:108-109
###################################################
warp.lsmobj@linfct


###################################################
### code chunk number 14: lsmeans-changes.rnw:118-127
###################################################
library(lme4)
data(Oats, package = "nlme")
Oats.lmer <- lmer(yield ~ factor(nitro) + Variety + (1|Block/Variety), 
    data = Oats, subset = -c(1,2,3,5,8,13,21,34,55))

### OLD
.old.lsmeans(Oats.lmer, pairwise ~ Variety)
### NEW
lsmeans(Oats.lmer, pairwise ~ Variety)


###################################################
### code chunk number 15: lsmeans-changes.rnw:133-137
###################################################
### OLD
.old.lsmeans(Oats.lmer, ~ nitro, at = list(nitro = c(.1,.2,.3)))
### NEW
lsmeans(Oats.lmer, ~ nitro, at = list(nitro = c(.1,.2,.3)))


###################################################
### code chunk number 16: lsmeans-changes.rnw:142-146
###################################################
(Oats.rg <- ref.grid(Oats.lmer))
Oats.quad <- update(Oats.lmer, yield ~ Variety + poly(nitro,2) + (1|Block/Variety))
ref.grid(Oats.quad)
ref.grid(Oats.quad, at = list(nitro = c(.1,.2,.3)))


###################################################
### code chunk number 17: lsmeans-changes.rnw:151-152
###################################################
(Oats.lsm <- lsmeans(Oats.rg, "nitro", by = "Variety"))


###################################################
### code chunk number 18: lsmeans-changes.rnw:155-157
###################################################
str(Oats.lsm)
(Oats.n <- lsmeans(Oats.lsm, "nitro"))


###################################################
### code chunk number 19: lsmeans-changes.rnw:162-163
###################################################
slotNames(Oats.lsm)


###################################################
### code chunk number 20: lsmeans-changes.rnw:169-170 (eval = FALSE)
###################################################
## summary(Oats.lsm, by = "nitro")


###################################################
### code chunk number 21: lsmeans-changes.rnw:173-174
###################################################
summary(Oats.n, infer = c(TRUE,TRUE))


###################################################
### code chunk number 22: lsmeans-changes.rnw:180-181
###################################################
(warp.con <- contrast(warp.lsmobj, method = "poly"))


###################################################
### code chunk number 23: lsmeans-changes.rnw:184-185
###################################################
contrast(warp.con, "revpairwise", by = "contrast")


###################################################
### code chunk number 24: lsmeans-changes.rnw:189-190
###################################################
cld(Oats.n, sort = FALSE)


###################################################
### code chunk number 25: lsmeans-changes.rnw:194-197
###################################################
chick.lmer <- lmer(weight ~ Time * Diet + (0 + Time | Chick), data = ChickWeight)
chick.lst <- lstrends(chick.lmer, ~ Diet, var = "Time")
cld(chick.lst, Letters = "|||||")


###################################################
### code chunk number 26: lsmeans-changes.rnw:206-207
###################################################
head(MOats)


###################################################
### code chunk number 27: lsmeans-changes.rnw:210-212
###################################################
MOats.mlm <- lm(yield ~ Block + Variety, data = MOats)
(MOats.rg <- ref.grid(MOats.mlm, mult.levs = list(nitro = c(0,.2,.4,.6))))


###################################################
### code chunk number 28: lsmeans-changes.rnw:217-219
###################################################
lsmeans(MOats.rg, ~ nitro)
lsmeans(MOats.rg, ~ Variety)


###################################################
### code chunk number 29: lsmeans-changes.rnw:222-224
###################################################
MOats <- transform(MOats, avg.yield = apply(yield, 1, mean))
lsmeans(lm(avg.yield ~ Block + Variety, data = MOats), ~ Variety)


