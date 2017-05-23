# functions that perform aspects of the MH regime
 library(MonoPoly)
 library(mvtnorm)

 source("basisGenerator.R")

# not needed after using BMonPol, although might still work!
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

  #sigma_mat <- diag(length(gamma_old)) * innov_sigma
  is_it_monotone_yet <- FALSE
  while (!is_it_monotone_yet) {
    gamma_proposal <- rnorm(length(gamma_old), mean = gamma_old, sd = innov_sigma)
    temp_beta <- gammaToBeta(t(t(gamma_proposal)), R_inv)
    # print(temp_beta
    # Hard coded to be over [0,1] at the moment, as we have rescaled data.
    is_it_monotone_yet <- ismonotone(temp_beta, a = 0, b = 1)
  }
  return(gamma_proposal)
}

genVarProp <- function(var_prev, innov_sd_var) {
  log_var_prev <- log(var_prev)
  log_var_prop <- rnorm(1, mean = log_var_prev, sd = innov_sd_var)
  return(exp(log_var_prop))
}


mvnPDF <-

calcAcceptProb <- function(Q_old, Q_prop, gamma_prop, gamma_old, var_prop,
                                var_old, y_vec, q_dim_swap) {
   # more conveinient to calculate mu as Q * gamma rather than X * Beta
 mu_new <- Q_prop %*% gamma_prop
 mu_old <- Q_old %*% gamma_old

   # comes from having flat priors on beta and sigma sq,
   # i.e. f(beta, sigma%2) \propto 1/sigma^2
 prior_ratio <- (1 / var_prop) / (1 / var_old)

 ratio <- dmvnorm(x = t(y_vec), mean = t(mu_new), sigma = var_prop * diag(length(y_vec))) /
          dmvnorm(x = t(y_vec), mean = t(mu_old), sigma = var_old * diag(length(y_vec)))

 return(ratio * prior_ratio * (1 / q_dim_swap))

}

