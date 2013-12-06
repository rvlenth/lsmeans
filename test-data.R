# Test cases for new lsmeans design

# subset of warpbreaks
ws = seq(2,53,by=3)

# add a fake covariate
mywarp = transform(warpbreaks, x = rnorm(54))

warp.data = lm(breaks ~ poly(x,3) + wool*tension, data = mywarp, subset = ws)
warp.log = lm(log(breaks) ~ poly(x,3) + wool*tension, data = mywarp, subset=ws)
warp.with = with(mywarp, lm(breaks ~ poly(x,3) + wool*tension, subset = ws))
attach(mywarp)
warp.att = lm(breaks ~ poly(x,3) + wool*tension, subset = ws)
detach()
