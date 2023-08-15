#' @export
tls.ROF <- function(Y, X, noise1, noise2, nsig, nsim.CI = 1000, nsim.rcc = 1000, REG = TRUE, df2 = NULL, plev = .9, rcc.flg = TRUE) {
  # Mandatory input: Y -- n x 1 matrix
  #        X -- n x nx matrix, nx is number of signals
  #        noise1 -- nnoise1 x n matrix, each row of this matrix is one piece of noise
  #        noise2 -- nnoise2 x n matrix, each row of this matrix is one piece of noise
  #        nsig -- vector, length of nx, ensemble size for calculating each signal
  # Optinal input:
  #        nsim.CI -- simulation size for confidence interval, default as 1000
  #        nsim.rcc -- simulation size for null distirbution of residual check stat
  #        df2  --  degree of freedom for noise2, if no input will treat noise2 pieces as
  #                 independent, namely, df2=nnoise2
  #        REG  -- regularization flag, apply regularization on noise1 or not
  #        rcc.flg -- flag to simulate empirical null distribution for residual test stat, default as TRUE
  # sample command lines:
  # x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
  # o1<-tls.ROF(x@Y,x@X,x@noise1,x@noise2,df2=180,REG=TRUE,nsig=5,rcc.flg=T)

  n <- length(Y)
  nx <- ncol(X)
  nn1 <- nrow(noise1)
  nn2 <- nrow(noise2)
  if (is.null(df2)) df2 <- nn2
  if (length(nsig) != nx) stop("tls.ROF: input nsig length not match input X column number")
  sqrt.matrix <- function(a) { # calculate square root of square matrix
    if (nrow(a) != ncol(a)) stop("sqrt.matrix: input matrix is not square")
    a.eig <- eigen(a)
    a.sqrt <- a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    return(a.sqrt)
  }

  in.sqrt.matrix <- function(a) { # calculate inverse square root of square matrix
    if (nrow(a) != ncol(a)) stop("in.sqrt.matrix: input matrix is not square")
    a.eig <- eigen(a)
    a.sqrt <- a.eig$vectors %*% diag(1 / sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    return(a.sqrt)
  }

  C1 <- if (REG) Creg(noise1) else t(noise1) %*% noise1 / nn1
  C12 <- t(in.sqrt.matrix(C1))
  gettlsbeta <- function(X, Y, C12, nsig) {
    nx <- ncol(X)
    Z <- cbind(C12 %*% (X %*% diag(sqrt(nsig), ncol = nx, nrow = nx)), C12 %*% Y)
    u <- svd(Z)
    ns <- ncol(X)
    nd <- length(u$d)
    v <- u$v[, nd]
    beta0 <- rep(NA, ns)
    for (i in 1:ns) beta0[i] <- -v[i] * sqrt(nsig[i]) / v[ns + 1]
    oout <- list()
    oout$u <- u
    oout$beta <- as.matrix(beta0, nrow = nx, ncol = 1)
    return(oout)
  }
  o0 <- gettlsbeta(X, Y, C12, nsig)
  beta0 <- o0$beta # finished estimate beta_hat
  u <- o0$u
  nd <- length(u$d)
  d <- u$d^2
  v <- t(u$v[, nd])
  dhat <- d
  n2t <- t(C12 %*% t(noise2))
  ui <- t(u$u[, nd])
  r1.stat <- d[nd] / ((ui %*% t(n2t) %*% n2t %*% t(ui)) / nn2) # residual test from A03
  r2.stat <- d[nd] / (sum(ui^2 %*% t(n2t^2)) / nn2) # residual test from ROF

  # now start estimate CI of beta, multi-signal cases need simulations for nsig sphere
  for (i in 1:(nx + 1)) {
    vi <- t(u$u[, i])
    dhat[i] <- d[i] / (vi %*% t(n2t) %*% n2t %*% t(vi) / nn2)
  }
  delta.dhat <- dhat - min(dhat)
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
  betaup <- betalow <- rep(NA, nx)
  vc.pts <- bi %*% t(u$v)
  for (i in 1:nx) vc.pts[, i] <- vc.pts[, i] * sqrt(nsig[i])
  for (j in 1:nx) {
    vc.d.pts <- complex(real = vc.pts[, j], imaginary = vc.pts[, (nx + 1)])
    vc.d.ref <- complex(real = v[j], imaginary = v[(nx + 1)])
    vprod.d <- vc.d.pts / vc.d.ref
    arg <- sort(Im(log(vprod.d)))
    delta.max1 <- max(arg[-1] - arg[-length(arg)])
    delta.max <- max(delta.max1, arg[1] - arg[length(arg)] + 2 * pi)
    arg.ref <- Im(log(vc.d.ref))
    arg.min <- min(arg) + arg.ref
    arg.max <- max(arg) + arg.ref
    betalow[j] <- -1 / tan(arg.min)
    betaup[j] <- -1 / tan(arg.max)
  }

  Tresi.simu.TLS <- function(Y, X, C, nsim = 1000, nn1, nn2, nsig) {
    n <- length(Y)
    nx <- ncol(X)
    C12 <- t(sqrt.matrix(C))
    C12.in <- t(in.sqrt.matrix(C))
    Z <- cbind(C12.in %*% (X %*% diag(sqrt(nsig), ncol = nx, nrow = nx)), C12.in %*% Y)
    u <- svd(Z)
    ns <- ncol(X)
    nd <- length(u$d)
    v <- u$v[, dim(u$v)[2]]
    beta0 <- as.matrix(-v[1:ns] * sqrt(nsig[1:ns]) / v[ns + 1], ncol = 1, nrow = nx)
    r1.stat <- r2.stat <- rep(NA, nsim)
    for (ith in c(1:nsim)) {
      Ys <- X %*% beta0 + C12 %*% rnorm(n, 0, 1)
      Xs <- X + C12 %*% (matrix(rnorm(n * nx), ncol = nx, byrow = T) %*% diag(1 / sqrt(nsig), ncol = nx, nrow = nx))
      noise1 <- t(C12 %*% matrix(rnorm(n * nn1, 0, 1), nrow = n, byrow = T))
      noise2 <- t(C12 %*% matrix(rnorm(n * nn2, 0, 1), nrow = n, byrow = T))
      C1hat <- Creg(noise1)
      Cs12 <- sqrt.matrix(C1hat)
      Zs <- cbind(C12.in %*% (Xs %*% diag(sqrt(nsig), ncol = nx, nrow = nx)), C12.in %*% Ys)
      u1 <- svd(Zs)
      nd <- length(u1$d)
      ids <- order(u1$d, decreasing = T)
      ui1 <- t(u1$u[, ids[nd]])
      v1 <- t(u1$v[, ids[nd]])
      n2t <- t(C12.in %*% t(noise2))
      d2 <- u1$d[ids[nd]]^2
      r1.stat[ith] <- d2 / ((ui1 %*% t(n2t) %*% n2t %*% t(ui1)) / nn2)
      r2.stat[ith] <- d2 / (sum(ui1^2 %*% t(n2t^2)) / nn2)
    }
    oout <- list()
    oout$r1.stat <- r1.stat
    oout$r2.stat <- r2.stat
    return(oout)
  }

  if (rcc.flg) {
    noisea <- rbind(noise1, noise2)
    Ca <- t(noisea) %*% noisea / (nn1 + nn2)
    rcc <- Tresi.simu.TLS(Y, X, Ca, nsim.rcc, nn1, nn2, nsig)
  }

  oout <- list()
  oout$beta <- beta0
  oout$r1.stat <- r1.stat
  oout$r2.stat <- r2.stat
  oout$n <- n
  oout$beta.CI <- cbind(betalow, beta0, betaup)
  if (rcc.flg) oout$rcc <- rcc
  return(oout)
}
