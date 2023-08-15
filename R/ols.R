#' @export
ols <- function(Y, X, noise1, noise2, nsig, df2 = NULL, plev = .95) {
  # sample command lines:
  # x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
  # o1<-ols(x@Y,x@X,x@noise1,x@noise2,df2=NULL)
  n <- length(Y)
  nx <- ncol(X)
  nn1 <- nrow(noise1)
  nn2 <- nrow(noise2)
  if (length(nsig) != nx) stop("ols: input nsig length not match input X column number")
  nmin <- min(c(n, nn1, nn2))
  if (is.null(df2)) df2 <- nn2
  if (df2 - nmin + 1 < 0) stop("df2 is to small to carry out residual check")
  C1 <- t(noise1) %*% noise1 / nn1
  n1eigen <- eigen(C1)
  ev1 <- n1eigen$values
  eof1 <- n1eigen$vector
  p <- eof1
  for (i in 1:ncol(p)) p[, i] <- eof1[, i] / max(sqrt(ev1[i]), .00005)
  # calculate second covariance matrix (cn2) from n2
  C2 <- t(noise2) %*% noise2 / nn2
  PX <- t(eof1) %*% X
  PY <- t(eof1) %*% Y
  noise2t <- noise2 %*% eof1
  output <- NULL
  for (ith in (nx + 1):nmin) {
    f <- solve(t(X) %*% p[, 1:ith] %*% t(p[, 1:ith]) %*% X) %*% t(X) %*% p[, 1:ith] %*% t(p[, 1:ith])
    tbeta <- f %*% Y # value of beta_hat
    vbeta <- f %*% C2 %*% t(f) # variance of beta_hat
    betalow <- tbeta - qt(plev, df2) * sqrt(diag(vbeta)) * sqrt((nsig + 1) / nsig)
    betaup <- tbeta + qt(plev, df2) * sqrt(diag(vbeta)) * sqrt((nsig + 1) / nsig)
    cn2t <- t(noise2t[, 1:ith]) %*% noise2t[, 1:ith] / nn2
    n2teigen <- eigen(cn2t)
    ev2 <- n2teigen$values
    ev2[ev2 < 0.00005] <- 0.00005
    eof2 <- n2teigen$vector
    u <- PY - PX %*% tbeta # residual base on ith truncation
    tmp <- rep(NA, ith)
    for (j in 1:ith) tmp[j] <- sum(u[1:ith] * eof2[1:ith, j])
    stat <- sum(tmp^2 / ev2[1:ith])
    dfn <- ith - nx
    # dfd, and modeifed F-stat proposed in Ribes13
    dfd <- nn2 - ith + 1
    f1 <- dfn * df2 * qf(.05, dfn, dfd) / (df2 - ith + 1)
    f2 <- dfn * df2 * qf(.95, dfn, dfd) / (df2 - ith + 1)
    output <- rbind(output, c(ith, as.vector(t(cbind(betalow, tbeta, betaup))), stat, f1, f2))
  }
  cnames <- "#EOF"
  for (i in 1:nx) cnames <- c(cnames, paste("beta", i, "_", c("low", "hat", "up"), sep = ""))
  cnames <- c(cnames, "RCstat", "RClow", "RCup")
  colnames(output) <- cnames
  return(output)
}
