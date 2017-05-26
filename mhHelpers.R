# functions that perform aspects of the MH regime
 library(MonoPoly)
 library(mvtnorm)

 source("basisGenerator.R")

# # not needed after using BMonPol, although might still work!
# genInitalGamma <- function(dat_df, d, R_inv, Q) {
#
#   Q_temp <- Q
#   R_inv_temp <- R_inv
#
#   # find the maximum uncostrained monotonic fit and then add zeros  to it
#   ii <- 0
#   min_model_found <- FALSE
#   while (!min_model_found) {
#
#     Q_temp <- Q_temp[, -(d  + 1 - ii - 1)]
#     R_inv_temp <- R_inv_temp[-(d  + 1 - ii - 1), -(d  + 1 - ii - 1)]
#
#     # dat_df$temp_fitted <- Q_temp %*% crossprod(Q_temp, dat_df$y)
#     # ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
#     #   geom_line(mapping = aes(x = x, y = temp_fitted))
#     #
#
#
#     pseudo_beta <- R_inv_temp %*% crossprod(Q_temp, dat_df$y)
#
#     min_model_found <- ismonotone(pseudo_beta, a = 0, b = 1)
#     ii <- ii + 1
#   }
#
#   gamma_proposal <- c(betaToGamma(pseudo_beta, R_inv_temp), rep(0, ii))
#
#   # checker <- FALSE
#   # while (!checker) {
#   #   gamma.proposal <- t(rmvnorm(n = 1, mean = gamma, sigma = diag(nrow = d + 1) / 1000 ))
#   #   beta.prop <- gammaToBeta(gamma.proposal, R_inv)
#   #   checker <- ismonotone(beta.prop)
#   # }
#   return(gamma_proposal)
# }


genMonoProp <- function(gamma_old, innov_sigma, R_inv) {
  # generate a random walk gamma proposal that is monotonic.

  #sigma_mat <- diag(length(gamma_old)) * innov_sigma
  is_it_monotone_yet <- FALSE
  while (!is_it_monotone_yet) {
    gamma_proposal <- rnorm(length(gamma_old), mean = gamma_old, sd = innov_sigma)
    temp_beta <- gammaToBeta(t(t(gamma_proposal)), R_inv)
    # Hard coded to be over [0,1] at the moment, as we have rescaled data.
    is_it_monotone_yet <- ismonotone(temp_beta, a = 0, b = 1)

  }

  # print("Gamma proposal is:")
  # print(gamma_proposal)

  # flog.debug("Gamma Proposal is: \n %f", gamma_proposal)

  return(gamma_proposal)
}

genVarProp <- function(var_prev, innov_sd_var) {
  log_var_prev <- log(var_prev)
  log_var_prop <- rnorm(1, mean = log_var_prev, sd = innov_sd_var)

  # flog.trace("Variance Proposal is %d", )

  # print("Variance proposal is:")
  # print(exp(log_var_prop))


  return(exp(log_var_prop))
}



calcPostRatio <- function(gamma_prop, gamma_curr, var_prop, var_curr, Q_full,
                          d_prop, d_curr, y_vec) {
  mu_prop <- Q_full[, 1:(d_prop + 1)] %*% gamma_prop
  mu_curr <- Q_full[, 1:(d_curr + 1)] %*% gamma_curr


  # comes from having flat priors on beta and sigma sq,
  #    # i.e. f(beta, sigma%2) \propto 1/sigma^2
  prior_ratio <- (1 / var_prop) / (1 / var_curr)



  ratio <- dmvnorm(x = t(y_vec), mean = t(mu_prop),
                   sigma = var_prop * diag(length(y_vec))) /
          dmvnorm(x = t(y_vec), mean = t(mu_curr),
                  sigma = var_curr * diag(length(y_vec)))

  # print("Post ratio is:")
  # print(ratio * prior_ratio)


  return(ratio * prior_ratio)
}

calcPropRatio <- function(gamma_k_prime_prop, k_prime_init, innov_sigma,
                          d_prop, d_curr, d_min, d_max) {

  if (d_prop == d_curr) {
    # if we don't propose a dimension change, then we don't have the dimension
    # jumping factor, nor the generation associated with the extra coefficient
    # so the proposal ratio is just 1 (everything else is random walk)
    return(1)
  }

  q_gamma_k_prime <- dnorm(x = gamma_k_prime_prop, mean = k_prime_init,
                           sd = innov_sigma)

  # print("werid Q value")
  # print(q_gamma_k_prime)

  dim_prop_ratio <- 1

  if ((d_prop == d_max) & (d_curr == (d_max - 1))) {
    there <- 1/3
    back <- 1/2
    dim_prop_ratio <- back / there

  } else if ((d_prop == (d_max - 1)) & (d_curr == d_max )) {
    there <- 1/2
    back <- 1/3
    dim_prop_ratio <- back / there

  } else if ((d_prop == d_min) & (d_curr == (d_min  + 1))) {
    there <- 1/3
    back <- 1/2
    dim_prop_ratio <- back / there

  } else if ((d_prop == (d_min + 1)) & (d_curr == d_min)) {
    there <- 1/2
    back <- 1/3
    dim_prop_ratio <- back / there

  }

  res <- (1 / (q_gamma_k_prime ^ 2)) * dim_prop_ratio

  # print("Proposal Ratio is:")
  # print(res)


  return(res)

}
