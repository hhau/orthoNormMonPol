# functions that perform aspects of the MH regime
 library(MonoPoly)
 library(mvtnorm)

 source("basisGenerator.R")

genInitalGamma <- function(dat_df, d, R_inv, Q) {

  Q_temp <- Q
  R_inv_temp <- R_inv

  # find the maximum uncostrained monotonic fit and then add zeros  to it
  ii <- 0
  min_model_found <- FALSE
  while (!min_model_found) {

    Q_temp <- Q_temp[, -(d  + 1 - ii - 1)]
    R_inv_temp <- R_inv_temp[-(d  + 1 - ii - 1), -(d  + 1 - ii - 1)]

    # dat_df$temp_fitted <- Q_temp %*% crossprod(Q_temp, dat_df$y)
    # ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
    #   geom_line(mapping = aes(x = x, y = temp_fitted))
    #


    pseudo_beta <- R_inv_temp %*% crossprod(Q_temp, dat_df$y)

    min_model_found <- ismonotone(pseudo_beta, a = 0, b = 1)
    ii <- ii + 1
  }

  gamma_proposal <- c(betaToGamma(pseudo_beta, R_inv_temp), rep(0, ii))

  # checker <- FALSE
  # while (!checker) {
  #   gamma.proposal <- t(rmvnorm(n = 1, mean = gamma, sigma = diag(nrow = d + 1) / 1000 ))
  #   beta.prop <- gammaToBeta(gamma.proposal, R_inv)
  #   checker <- ismonotone(beta.prop)
  # }
  return(gamma_proposal)
}


genMonoProp <- function(gamma_old, innov_sigma, R_inv) {
  # generate a random walk gamma proposal that is monotonic.

  sigma_mat <- diag(length(gamma_old)) * innov_sigma
  is_it_monotone_yet <- FALSE
  while (!is_it_monotone_yet) {
    gamma_proposal <- rmvnorm(1, mean = gamma_old, sigma = sigma_mat)
    temp_beta <- gammaToBeta(t(gamma_proposal), R_inv)
    print(temp_beta)
    is_it_monotone_yet <- ismonotone(temp_beta)
  }
  return(gamma_proposal)
}

calcAcceptProb <- function() {

}
