# Mandatory input: Y -- n x 1 matrix
#        X -- n x nx matrix, nx is number of signals
#        noise1 -- nnoise1 x n matrix, each row of this matrix is one piece of noise
#        noise2 -- nnoise2 x n matrix, each row of this matrix is one piece of noise
#        nsig -- vector, length of nx, ensemble size for each signal
# Optional input:
#        nsim.CI -- simulation size for confidence interval, default as 1000
#        df2  --  degree of freedom for noise2, if no input will treat noise2 pieces as
#                 independent, namely, df2=nnoise2
#        REG  -- regularization flag, apply regularization on noise1 or not, default as FALSE
#        plev -- confidence level used in residual check, defaul as 0.9
# output as matrix (n-nx) rows (corresponding to EOF number), each row contains EOF#, (beta_low,
#        beta_hat, beta_up) for every signal
# sample command lines:
# x <- readin('obs_sig.txt','noise1.txt','noise2.txt')
# o1 <- tls.A03(x@Y,x@X,x@noise1,x@noise2,nsig=5,df2=220,REG=TRUE)

#' @export 
tls.A03 <- function(Y, X, noise1, noise2, nsig, nsim.CI = 1000, df2 = NULL, REG = FALSE, plev = .9) {
  
  if (missing(Y)) stop("input Y missing in tls.A03")
  if (missing(X)) stop("input X missing in tls.A03")
  if (missing(noise1)) stop("input noise1 missing in tls.A03")
  if (missing(noise2)) stop("input noise2 missing in tls.A03")
  if (missing(nsig)) stop("input nsig missing in tls.A03")

  n <- length(Y)
  nx <- ncol(X)
  nn1 <- nrow(noise1)
  nn2 <- nrow(noise2)
  if (is.null(df2)) df2 <- nn2
  if (length(nsig) != nx) stop("tls.A03: input nsig length not match input X column number")
  if (REG) C1 <- Creg(noise1) else C1 <- t(noise1) %*% noise1 / nn1
  x <- X
  for (i in 1:nx) x[, i] <- X[, i] * sqrt(nsig[i])
  # x is adjusted X with nsig
  eigen1 <- eigen(C1)
  P <- eigen1$vectors
  for (i in 1:min(n, nn1)) P[, i] <- P[, i] / sqrt(eigen1$values[i])
  P <- t(P)
  # P is pre-whitening operator
  Z <- cbind(P %*% x, P %*% Y)
  betaout <- NULL
  betaCI <- NULL
  for (ith in (nx + 1):n) {
    zz <- Z[1:ith, ]
    u <- svd(zz)
    betahat <- -u$v[1:nx, (nx + 1)] / u$v[(nx + 1), (nx + 1)] * sqrt(nsig)
    d <- u$d^2
    nd <- length(u$d)
    ui <- t(u$u[, nd])
    n2t <- noise2 %*% t(P[1:ith, ])
    r1.stat <- d[nd] / ((ui %*% t(n2t) %*% n2t %*% t(ui)) / nn2)
    dhat <- rep(NA, nx + 1)
    for (i in 1:(nx + 1)) {
      vi <- t(u$u[, i])
      dhat[i] <- d[i] / (vi %*% t(n2t) %*% n2t %*% t(vi) / nn2)
    }
    delta.dhat <- dhat - min(dhat)
    # sampling on critical ellipse to construct correspoding beta_simus,
    # then get max/min value as upper/lower bound for scaling factor beta_hat
    Crit <- qf(plev, 1, df2)
    if (nx > 1) {
      unit.tmp <- matrix(rnorm(nsim.CI * nx, 0, 1), nrow = nsim.CI, byrow = T)
      unit <- t(apply(unit.tmp, 1, function(x) {
        x / sqrt(sum(x^2))
      }))
    } else {
      unit <- t(t(c(1, -1)))
    }
    ai <- unit * sqrt(Crit)
    bi <- cbind(ai, NA)
    for (i in 1:nx) bi[, i] <- ai[, i] / sqrt(delta.dhat[i])
    bi[, nx + 1] <- apply(t(t(bi[, 1:nx])), 1, function(x) {
      sqrt(1 - sum(x^2))
    })
    nsim <- dim(bi)[1]
    betaup <- betalow <- betaup.new <- betalow.new <- rep(NA, nx)
    vc.pts <- bi %*% t(u$v)
    for (i in 1:nx) vc.pts[, i] <- vc.pts[, i] * sqrt(nsig[i])
    for (j in 1:nx) {
      beta.tmp <- -vc.pts[, j] / vc.pts[, (nx + 1)]
      betaup[j] <- max(beta.tmp, na.rm = T)
      betalow[j] <- min(beta.tmp, na.rm = T)
    }
    obeta <- NULL
    for (i in 1:nx) obeta <- c(obeta, betalow[i], betahat[i], betaup[i])
    betaout <- rbind(betaout, c(ith, obeta, r1.stat, (ith - nx) * qf(.05, ith - nx, df2), (ith - nx) * qf(.95, ith - nx, df2)))
  }
  cnames <- "#EOF"
  for (i in 1:nx) cnames <- c(cnames, paste("beta", i, "_", c("low", "hat", "up"), sep = ""))
  cnames <- c(cnames, "RCstat", "RClow", "RCup")
  colnames(betaout) <- cnames
  return(betaout)
}
