plotbetas <- function(betaout, ...) {
  cnams <- colnames(betaout)
  sid <- grepl("_hat", cnams)
  sigs <- cnams[sid]
  nsig <- length(sigs)
  sigs.low <- gsub("hat", "low", sigs)
  sigs.up <- gsub("hat", "up", sigs)
  for (ith in 1:nsig) {
    yrange <- range(betaout[, c(sigs.low[ith], sigs.up[ith])], na.rm = T)
    yrange[1] <- min(0, yrange[1])
    yrange[2] <- max(yrange[2], 1)
    plot(betaout[, "#EOF"], betaout[, sigs[ith]],
      type = "p", main = paste("best estimates of scaling factors for ",
        strsplit(sigs[ith], "_")[[1]][1],
        sep = ""
      ), ylim = yrange,
      xlab = "Number of EOF patterns retained in the truncation", ylab = "scaling factors", col = "red", ...
    )
    abline(h = 0)
    abline(h = 1, lty = 3)
    for (j in 1:nrow(betaout)) {
      lines(rep(betaout[j, "#EOF"], 2), c(betaout[j, sigs.low[ith]], betaout[j, sigs.up[ith]]))
      lines(betaout[j, "#EOF"] + c(-.1, .1), rep(betaout[j, sigs.low[ith]], 2))
      lines(betaout[j, "#EOF"] + c(-.1, .1), rep(betaout[j, sigs.up[ith]], 2))
    }
  }
}

plotrstat <- function(betaout, ...) {
  yrange <- range(1 / betaout[, c("RCstat", "RClow", "RCup")], na.rm = T)
  plot(betaout[, "#EOF"], 1 / betaout[, "RCstat"],
    type = "p", col = "red", ylim = yrange,
    main = "residual consistency test",
    xlab = "Number of EOF patterns retained in the truncation",
    ylab = "Cumulative ratio model/observation variance", log = "y"
  )
  lines(betaout[, "#EOF"], 1 / betaout[, "RClow"], type = "l", lty = 3, col = "blue")
  lines(betaout[, "#EOF"], 1 / betaout[, "RCup"], type = "l", lty = 3, col = "blue")
}
