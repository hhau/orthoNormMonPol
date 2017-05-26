# functions specifically to enable the Reversible jump aspect of this routine
source("mhHelpers.R")

dimProposer <- function(d_current, d_min, d_max) {
  # currently all these have symmetric proposals
  # so as to not increase the number of terms in my acceptance ratio
  # i.e.
  # q(k' - > k) = q(k - > k')

  if (d_current < d_min | d_current > d_max) {
   stop("current dimension is outside of allowed bounds")
  }

  if (d_current == d_min) {
    res <- sample(x = c(d_current, d_current + 1), size = 1)

  } else if (d_current == d_max) {
    res <- sample(x = c(d_current, d_current - 1), size = 1)

  } else {
    res <- sample(x = c(d_current - 1, d_current, d_current + 1), size = 1)

  }

  return(res)
}

genGammaProp <- function(current_gamma, d_prop, d_curr,
                         innov_sigma, R_inv_full, full_inits) {
  # make sure the full R_inv gets passed to this function, so we can pass
  # the correct dimension down

  if (d_prop == d_curr) {
    gamma_prop <- genMonoProp(gamma_old = current_gamma,
                              innov_sigma = innov_sigma,
                              R_inv = R_inv_full[1:(d_curr + 1), 1:(d_curr + 1)])

  } else if (d_prop ==  (1 + d_curr)) {
    prop_mean <- c(current_gamma, full_inits[d_prop + 1])
    gamma_prop <- genMonoProp(gamma_old = prop_mean, innov_sigma = innov_sigma,
                              R_inv = R_inv_full[1:(d_prop + 1),1:(d_prop + 1)])

  } else if (d_prop == (d_curr - 1)) {
    prop_mean <- current_gamma[1:(d_prop + 1)]
    gamma_prop <- genMonoProp(gamma_old = prop_mean, innov_sigma = innov_sigma,
                              R_inv = R_inv_full[1:(d_prop + 1), 1:(d_prop + 1)])

  }

  # print("Double checking what gamma_proposal is:")
  # print(gamma_prop)


  return(gamma_prop)

}

calcAcceptProb <- function(gamma_prop, gamma_curr, d_prop, d_curr, d_min, d_max,
                           var_prop, var_curr, Q_full, full_inits, y_vec, innov_sigma) {

  accept_prob <- NULL

  if (d_prop == d_curr) {
    post <- calcPostRatio(gamma_prop = gamma_prop, gamma_curr = gamma_curr,
                          var_prop = var_prop, var_curr = var_curr,
                          Q_full = Q_full, d_prop = d_prop, d_curr = d_curr,
                          y_vec = y_vec)
    prop <- 1 # could put the function call here but why?

    accept_prob <- min(1, post * prop)

  } else if (d_prop == (d_curr + 1)) {

    post <- calcPostRatio(gamma_prop = gamma_prop, gamma_curr = gamma_curr,
                          var_prop = var_prop, var_curr = var_curr,
                          Q_full = Q_full, d_prop = d_prop, d_curr = d_curr,
                          y_vec = y_vec)

    prop <- calcPropRatio(gamma_k_prime_prop = gamma_prop[d_prop + 1],
                          k_prime_init = full_inits[d_prop + 1],
                          innov_sigma = innov_sigma, d_prop = d_prop, d_curr = d_curr,
                          d_min = d_min, d_max = d_max)
    accept_prob <- min(1, post * prop)

  } else if (d_prop == (d_curr - 1)) {
     # different here, have to ~flip it n reverse it~
     # basically have to pretend proposal is current, and current is proposal
     # then invert after calculating accept_prob

    post <- calcPostRatio(gamma_prop = gamma_curr,
                          gamma_curr = gamma_prop,
                          var_prop = var_curr,
                          var_curr = var_prop,
                          Q_full = Q_full,
                          d_prop = d_curr,
                          d_curr = d_prop,
                          y_vec = y_vec)

    prop <- calcPropRatio(gamma_k_prime_prop = gamma_curr[d_curr + 1],
                          k_prime_init = full_inits[d_curr + 1],
                          innov_sigma = innov_sigma,
                          d_prop = d_curr,
                          d_curr = d_prop,
                          d_min = d_min, d_max = d_max)
    temp <- (post * prop)^(-1)
    accept_prob <- min(1, temp)

  }

  # print("Acceptance probability is:")
  # print(accept_prob)

  return(accept_prob)

}
