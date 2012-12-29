### functions to implement different families of contrasts
### All return a matrix or data frame whose columns are the desired contrasts coefs
### with appropriate row and column names
### Also they have two attributes: 
###   "desc" is an expanded description of the family,
###   "adjust" is the default multiplicity adjustment (used if adjust="auto" in lsmeans)

# all pairwise trt[i] - trt[j], i < j
pairwise.lsmc = function(levs,...) {
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in 1:(k-1)) {
        for (j in (i+1):k) {
            con = rep(0,k)
            con[i] = 1
            con[j] = -1
            nm = paste(levs[i], levs[j], sep = " - ")
            M[[nm]] = con
        }
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "pairwise differences"
    attr(M, "adjust") = "tukey"
    M
}

# all pairwise trt[j] - trt[i], j > i
revpairwise.lsmc = function(levs,...) {
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in 2:k) {
        for (j in 1:(i-1)) {
            con = rep(0,k)
            con[i] = 1
            con[j] = -1
            nm = paste(levs[i], levs[j], sep = " - ")
            M[[nm]] = con
        }
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "pairwise differences"
    attr(M, "adjust") = "tukey"
    M
}

# Poly contrasts - scaled w/ integer levels like most tables
# ad hoc scaling works for up to 13 levels
poly.lsmc = function(levs, max.degree=min(6,k-1)) {
    nm = c("linear", "quadratic", "cubic", "quartic", paste("degree",5:20))
    k = length(levs)
    M = as.data.frame(poly(1:k, min(20,max.degree)))
    for (j in 1:ncol(M)) {
        con = M[ ,j]
        pos = which(con > .01)
        con = con / min(con[pos])
        z = max(abs(con - round(con)))
        while (z > .05) {
            con = con / z
            z = max(abs(con - round(con)))
        }
        M[ ,j] = round(con)
    }
    row.names(M) = levs
    names(M) = nm[1:ncol(M)]
    attr(M, "desc") = "polynomial contrasts"
    attr(M, "adjust") = "none"
    M
}

# All comparisons with a control; ref = index of control group
# New version -- allows more than one control group (ref is a vector)
trt.vs.ctrl.lsmc = function(levs, ref=1) {
    if ((min(ref) < 1) || (max(ref) > length(levs)))
        stop("Reference levels are out of range")
    k = length(levs)
    cnm = ifelse(length(ref)==1, 
        levs[ref], 
        paste("avg(", paste(levs[ref], collapse=","), ")", sep=""))
    templ = rep(0, length(levs))
    templ[ref] = -1 / length(ref)
    M = data.frame(levs=levs)
    for (i in 1:k) {
        if (i %in% ref) next
        con = templ
        con[i] = 1
        nm = paste(levs[i], cnm, sep = " - ")
        M[[nm]] = con
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "differences from control"
    attr(M, "adjust") = "sidak"
    M
}

# control is 1st level
trt.vs.ctrl1.lsmc = function(levs, ...) {
    trt.vs.ctrl.lsmc(levs, ref=1)
}

# control is last level
trt.vs.ctrlk.lsmc = function(levs, ...) {
    trt.vs.ctrl.lsmc(levs, ref=length(levs))
}

