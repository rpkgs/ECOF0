setClass(
  Class = "ECOF",
  representation = representation(
    X = "matrix",
    Y = "matrix",
    noise1 = "matrix",
    noise2 = "matrix"
  )
)

#' @export
print.ECOF <- function(x, ...) str(x)

#' @export
setMethod("show", signature(object = "ECOF"), function(object) print(object))

# all the required data will be included in class tls, and corresponding
# validate function: checkOF() will apply checking rules when claim a
# new class. Later on, user might use different function: tls or tls_rof
# using tls class as input data, and also parameters for each function.
checkOF <- function(object) {
  if (length(object@X) == 0) {
    return("input X is empty")
  }
  if (length(object@Y) == 0) {
    return("input Y is empty")
  }
  if (length(object@Y) == 0) {
    return("input Y is empty")
  }
  if (length(object@Y) == 0) {
    return("input Y is empty")
  }
  if (ncol(object@Y) != 1) {
    return("Y should contains 1 column only")
  }
  np <- nrow(object@X)
  nx <- ncol(object@X)
  if (nrow(object@Y) != np) {
    return("input Y and X sample size differ")
  }
  if (ncol(object@noise1) != np) {
    return("input noise1 column number not equal observation number")
  }
  if (ncol(object@noise2) != np) {
    return("input noise2 column number not equal observation number")
  }
  return(TRUE)
}
setValidity("ECOF", checkOF)

#' @export
readin <- function(file.obssig, file.noise1, file.noise2) {
  # read from data files and return class of tls
  # sample command lines:
  # x<-readin('obs_sig.txt','noise1.txt','noise2.txt')
  # x will be input for ols, tls or tls_rof
  itmp <- try(obssig <- as.matrix(read.table(file.obssig)), T)
  if (inherits(itmp, "try-error")) stop(paste("read file:", file.obssig, "error"))
  Y <- t(t(obssig[1, ]))
  X <- if (nrow(obssig) > 2) t(obssig[-1, ]) else t(t(obssig[-1, ]))
  itmp <- try(noise1 <- as.matrix(read.table(file.noise1)), T)
  if (inherits(itmp, "try-error")) stop(paste("read file:", file.noise1, "error"))
  itmp <- try(noise2 <- as.matrix(read.table(file.noise2)), T)
  if (inherits(itmp, "try-error")) stop(paste("read file:", file.noise2, "error"))
  np <- length(Y)
  if (ncol(noise1) != np) stop("dim of noise1 not match obssig")
  if (ncol(noise2) != np) stop("dim of noise2 not match obssig")
  return(new(Class = "ECOF", X = X, Y = Y, noise1 = noise1, noise2 = noise2))
}
