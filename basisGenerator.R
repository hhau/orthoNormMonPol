# Supporting functions for OrthoNormalMonPol
# AAM / hhau 27/04/17

library(MonoPoly)


genAllMatrices <- function(x, d) {
  # generate the q, r_inverse matricies from a (column) vector of x values,
  # and a polynomial degree d.
  # Berwin's adaptation of ****** book notes
  x <- sort(x)

  # set up matrix w/ correct dimensions
  Q <- matrix(NA, nrow = length(x), ncol = d + 1)
  R_inv <- matrix(NA, nrow = d + 1, ncol = d + 1)

  # renormalising constants
  a <- b <- c <- 0
  gamma <- rep(0, d + 1)
  k <- rep(1, d + 1)

  # setup first two cols, corresponding to the 1 and x cols in the X matrix
  n <- length(unique(x))
  Q[ ,1] <- 1 / sqrt(n)
  R_inv[1, 1] <- k[1]

  Q[, 2] <- x - mean(x)
  k[2] <- 1 / (sqrt(sum(Q[, 2]^2)))
  Q[, 2] <- Q[, 2] * k[2]
  R_inv[1, 2] <- -mean(x)
  R_inv[2, 2] <- 1
  R_inv[, 2] <- R_inv[, 2] * k[2]

  gamma[1] <- k[2] / k[1]

  for (ii in seq_len(d - 1) + 2) {
    gamma[ii - 1] <- k[ii] / k[ii - 1]
    A <- -gamma[ii - 1] * sum(x*Q[, ii - 1]^2)
    B <- gamma[ii - 1]
    C <- (gamma[ii - 1] * sum(Q[, ii - 1]^2)) / (gamma[ii - 2] * sum(Q[, ii - 2]^2))

    Q[, ii] <- (A + B*x)*Q[,ii - 1] - C*Q[, ii - 2]
    ix1 <- 1:(ii - 1)
    ix <-  1:ii
    R_inv[ix, ii] <- A*R_inv[ix, ii - 1] + B*c(0, R_inv[ix1, ii - 1]) - C*R_inv[ix, ii - 2]

    ## reorthogonalisation
    for (jj in 1:(ii - 1)) {
      tmp <- crossprod(Q[, jj], Q[, ii])
      Q[, ii] <- Q[, ii] - tmp*Q[, jj]
      R_inv[,ii] <- R_inv[,ii] - tmp*R_inv[,jj]
    }

    k[ii] <- 1 / sqrt(sum(Q[, ii]^2))
    gamma[ii - 1] <- k[ii] / k[ii - 1]
    Q[, ii] <- Q[, ii] * k[ii]
    R_inv[, ii] <- R_inv[, ii] * k[ii]

  }
  return(list(Q = Q, R_inv = R_inv))
}
