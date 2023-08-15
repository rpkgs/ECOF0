# This ECOF package is coded by Feng, Yang (yang.feng@ec.gc.ca) with technical support and scientific review by Qiuzi Han Wen(HanQiuzi.Wen@ec.gc.ca).

# Main functions
# This package contains R functions to conduct detection and attribution analysis using three different algorithms under the optimal fingerprint framework. These include:
# 1.ols(), the Ordinary Least Squares method (Allen and Tett, 1999).
# 2.tls (), the Total Least Squares method (Allen and Scott, 2003). Note that the confidence intervals for the scaling factors are obtained using the method provided in the ROF package by Dr. A. Ribes (Ribes, A., 2012)
# 3.tls_rof (), the Regularized Optimal Fingerprint method (Ribes et al, 2013a). This function is translated from routines in the ROF package V0.8 coded by Dr. A. Ribes, using SCILAB. (Ribes, A., 2012)

Creg <- function(Cn) {
  # regularization for noise structure
  # input Cn is an n x p matrix, n is number of pieces of noise;
  # nx is point number for each piece
  nx <- ncol(Cn)
  n <- nrow(Cn)
  CE <- t(Cn) %*% Cn / n
  Ip <- diag(1, nx, nx)
  m <- mean(diag(CE))
  XP <- CE - diag(m, nx, nx)
  d2 <- mean(diag(XP %*% t(XP)))
  bt <- rep(NA, n)
  for (i in 1:n) {
    Mi <- t(t(Cn[i, ])) %*% t(Cn[i, ]) # Mi is nx x nx matrix
    bt[i] <- mean(diag((Mi - CE) %*% t(Mi - CE)))
  }
  bb2 <- sum(bt) / n^2
  b2 <- min(bb2, d2)
  a2 <- d2 - b2
  Creg <- b2 * m / d2 * Ip + a2 / d2 * CE
  return(Creg)
}

# given space-temporal structure and p1,p2 as base period start and end point
# location, return EOF(,nt-1) for reduce dimension operator

#' @export
redop <- function(nt, p1, p2) {
  if (p1 < 1 | p1 > nt) stop("input p1 error in reduce")
  if (p2 < 1 | p2 > nt) stop("input p2 error in reduce")
  M <- diag(1, nt)
  M[, p1:p2] <- M[, p1:p2] - 1 / (p2 - p1 + 1)
  eM <- eigen(M)
  if (abs(eM$values[nt]) > 1E-5) stop("ev error in reduce")
  if (any(abs(eM$values[-nt]) < 1E-5)) stop("ev error2 in reduce")
  u <- eM$vector[, -nt]
  return(u)
}

# given reduce dimension operator u and data vector x and space-temporal structure nt,ns
# return reduced vector (nt-1)*ns
redvec <- function(u, x, nt, ns, timefirst = T) {
  if (length(x) != nt * ns) stop("input dim error in redvec")
  if (all.equal(dim(u), c(nt, nt - 1)) != T) stop("input u dim error in redvec")
  x1 <- matrix(x, nrow = nt, byrow = timefirst)
  ox <- t(u) %*% x1
  return(as.vector(t(ox)))
}
# usage:
# for example, vector of obs length as 50, listed as 10yrs x 5points (data 1~5 for year 1,
#     6~10 for year 2, so on so for, this kind of sequence should apply to all X and noise
#     pieces), then nt=10, ns=5
#     if the data was anomanies wrt year 6~9, namely baseperiod as 6~9, then p1=6, p2=9
# u<-redop(nt=10,ns=5,p1=6,p2=9)
# newY<-apply(x@Y,2,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# newX<-apply(x@X,2,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# newnoise1<-apply(x@noise1,1,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)
# newnoise2<-apply(x@noise2,1,function(x,u=u,nt=10,ns=5){return(redvec(u,x,nt,ns))},u=u)

# redECOF is a function to reduce everything (Y,X,noise1,noise2) from a ECOF object and return
# a new ECOF object contains reduced elements (Y,X,noise1,noise2)
#' @export
redECOF <- function(x, u, nt, ns, timefirst = T) {
  newY <- apply(x@Y, 2, function(x, u, nt, ns, timefirst) {
    return(redvec(u = u, x, nt = nt, ns = ns, timefirst = timefirst))
  }, u = u, nt = nt, ns = ns, timefirst = timefirst)
  newX <- apply(x@X, 2, function(x, u, nt, ns, timefirst) {
    return(redvec(u = u, x, nt = nt, ns = ns, timefirst = timefirst))
  }, u = u, nt = nt, ns = ns, timefirst = timefirst)
  newnoise1 <- t(as.matrix(apply(x@noise1, 1, function(x, u, nt, ns, timefirst) {
    return(redvec(u = u, x, nt = nt, ns = ns, timefirst = timefirst))
  }, u = u, nt = nt, ns = ns, timefirst = timefirst)))
  newnoise2 <- t(as.matrix(apply(x@noise2, 1, function(x, u, nt, ns, timefirst) {
    return(redvec(u = u, x, nt = nt, ns = ns, timefirst = timefirst))
  }, u = u, nt = nt, ns = ns, timefirst = timefirst)))
  return(new(Class = "ECOF", X = newX, Y = newY, noise1 = newnoise1, noise2 = newnoise2))
}
