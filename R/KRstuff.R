### Kenward-Roger df for a 1-dimensional hypothesis
### Code borrowed and adapted from package pbkrtest.
### In particular, the function .KR_Adjust

# Non-exported utility fcns...
.spur = function(U){
  sum(diag(U))
}
.divZero = function(x,y,tol=1e-14){
  ## ratio x/y is set to 1 if both |x| and |y| are below tol
  x.y  =  if( abs(x)<tol & abs(y)<tol) {1} else x/y
  x.y
}

### Returns denom d.f. for testing lcoefs'beta = 0 where lcoefs is a vector
# PhiA is result of call to PhiA = vcovAdj(object, 0) in pbkrtest package
# Phi is vcov(object) ## May not now be needed
# lcoefs is contrast of interest
# varlb is my already-computed value of lcoef' Phi lcoef = est variance of lcoef'betahat
.KRdf.mer = function(PhiA, Phi, lcoefs, varlb) {
  # I guess I'm supposed to use unadjusted varlb, as much as I'm inclined not to
      vlb = sum(lcoefs * (Phi %*% lcoefs))
      Theta = Matrix(as.numeric(outer(lcoefs,lcoefs) / vlb), nrow=length(lcoefs))
      P = attr(PhiA,"P")
      W = attr(PhiA,"W")

      A1 = A2 = 0
      ThetaPhi = Theta%*%Phi
      n.ggamma = length(P)
      for (ii in 1:n.ggamma) {
        for (jj in c(ii:n.ggamma)) {
          e = ifelse(ii==jj, 1, 2)
          ui = ThetaPhi %*% P[[ii]] %*% Phi
          uj = ThetaPhi %*% P[[jj]] %*% Phi
          A1 =  A1 +  e* W[ii,jj] * (.spur(ui) * .spur(uj))
          A2 =  A2 +  e* W[ii,jj] *  sum(ui * t(uj))
        }}
    
# substituted q = 1 in pbkrtest code and simplified
      B  =  (A1 + 6 * A2) / 2
      g  =  (2 * A1 - 5 * A2)  / (3 * A2)
      c1 =  g/(3 + 2 * (1 - g))
      c2 =  (1 - g) / (3 + 2 * (1 - g))
      c3 =  (3 - g) / (3 + 2 * (1 - g))
      EE =  1 + A2
      VV =  2 * (1 + B)
      EEstar  =  1/(1 - A2)
      VVstar  =  2 * ((1 + c1 * B)/((1 - c2 * B)^2  *  (1 - c3 * B)))
      V0 = 1 + c1 * B
      V1 = 1 - c2 * B
      V2 = 1 - c3 * B
      V0 = ifelse(abs(V0) < 1e-10, 0, V0)
      rho  = (.divZero(1 - A2, V1))^2 * V0/V2
      df2  =  4 + 3 / (rho - 1)
      df2
}
