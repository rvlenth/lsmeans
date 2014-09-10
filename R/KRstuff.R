### Kenward-Roger df for a 1-dimensional hypothesis
### Code borrowed and adapted from package pbkrtest.
### In particular, the function .KR_adjust

### We will phase this out and use pbkrtest::Lb_ddf instead when it's available

# Non-exported utility fcns...
.spur = function(U){
  sum(diag(U))
}
.divZero = function(x,y,tol=1e-14){
  ## ratio x/y is set to 1 if both |x| and |y| are below tol
  x.y  =  if( abs(x)<tol & abs(y)<tol) {1} else x/y
  x.y
}


# New version that also accepts a matrix
# I just borrowed code from pbkrtest:::.KR_adjust, with slight changes
.KRdf.mer = function(PhiA, Phi, L) {
    if (!is.matrix(L))
        L = matrix(L, nrow = 1)
    Theta <- t(L) %*% solve(L %*% Phi %*% t(L), L)
    P <- attr(PhiA, "P")
    W <- attr(PhiA, "W")
    A1 <- A2 <- 0
    ThetaPhi <- Theta %*% Phi
    n.ggamma <- length(P)
    for (ii in 1:n.ggamma) {
        for (jj in c(ii:n.ggamma)) {
            e <- ifelse(ii == jj, 1, 2)
            ui <- ThetaPhi %*% P[[ii]] %*% Phi
            uj <- ThetaPhi %*% P[[jj]] %*% Phi
            A1 <- A1 + e * W[ii, jj] * (.spur(ui) * .spur(uj))
            A2 <- A2 + e * W[ii, jj] * sum(ui * t(uj))
        }
    }
    q <- nrow(L)        # instead of finding rank
    B <- (1/(2 * q)) * (A1 + 6 * A2)
    g <- ((q + 1) * A1 - (q + 4) * A2)/((q + 2) * A2)
    c1 <- g/(3 * q + 2 * (1 - g))
    c2 <- (q - g)/(3 * q + 2 * (1 - g))
    c3 <- (q + 2 - g)/(3 * q + 2 * (1 - g))
    EE <- 1 + (A2/q)
    VV <- (2/q) * (1 + B)
    EEstar <- 1/(1 - A2/q)
    VVstar <- (2/q) * ((1 + c1 * B)/((1 - c2 * B)^2 * (1 - c3 * B)))
    V0 <- 1 + c1 * B
    V1 <- 1 - c2 * B
    V2 <- 1 - c3 * B
    V0 <- ifelse(abs(V0) < 1e-10, 0, V0)
    rho <- 1/q * (.divZero(1 - A2/q, V1))^2 * V0/V2
    df2 <- 4 + (q + 2)/(q * rho - 1)
    df2
}


