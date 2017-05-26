# AAM / hhau 27/04/17
#
# Orhtonormal Bayesian Monotonic polynomial using metropolis hastings
# and a renormalised QR decomposition
#
# This is intending as a starting point for, upon which we will later add
# the ability to traverse polynomial degree via reversible jump algorithm.


# setup will probably look something like this
#
# load data /librarys / other files
#
# setup data structures for samples / inital values
#
# # then within the MH loop
# # generate a sample (need function for this), for the orthonormal parameters
# # theta
#
# # then accept reject that sample with is.monotonic after converting to
# # monomial basis
#
# # then calculate the MH acceptance proposal on that.

library(MonoPoly)
library(fda)
library(ggplot2)
library(mvtnorm)
library(rstan)
library(TSA)
library(futile.logger)



source("basisGenerator.R")
source("mhHelpers.R")
source("plotDiagnoseFunctions.R")
source("rjHelper.R")
source("dumb.R")

source("./BMonPolCode/getBMonPolBeta.R")


clearAndSetupSpace()

for (ii in 2:(n_iter + 1)) {

  d <- d_mat[ii - 1]
  d_prop <- dimProposer(d, d_min, d_max)

  # print("previous dim:")
  # print(d)
  #
  # print("proposed dim:")
  # print(d_prop)



  variance_proposal <- genVarProp(var_prev = var_mat[ii - 1],
                                  innov_sd_var = innov_sd_var)

  gamma_proposal <- genGammaProp(current_gamma = gamma_samples_list[[ii - 1]],
                                 d_prop = d_prop, d_curr = d, innov_sigma = innov_sd_beta,
                                 R_inv_full = R_inv_full, full_inits = gamma_inital)


  accept_probabililty <- calcAcceptProb(gamma_prop = gamma_proposal,
                                        gamma_curr = gamma_samples_list[[ii - 1]],
                                        d_prop = d_prop, d_curr = d, d_min = d_min,
                                        d_max = d_max, var_prop = variance_proposal,
                                        var_curr = var_mat[ii - 1], Q_full = Q_full,
                                        full_inits = gamma_inital, y_vec = y,
                                        innov_sigma = innov_sd_beta)






  if (runif(1) < accept_probabililty) {
    # gamma_mat[ii, ] <- gamma_proposal

    gamma_samples_list[[ii]] <- gamma_proposal

    d_mat[ii] <- d_prop


    var_mat[ii] <- variance_proposal
    n_accept <- n_accept + 1

    gamma_accepted_samples[[n_accept]] <- gamma_proposal
    d_accepted[[n_accept]] <- d_prop

    # print("Accepted")

  } else {
    # gamma_mat[ii, ] <- gamma_mat[ii - 1, ]
    var_mat[ii] <- var_mat[ii - 1]

    gamma_samples_list[[ii]] <- gamma_samples_list[[ii - 1]]

    d_mat[ii] <- d_mat[ii - 1]

  }

  if (ii %% 1000 == 0) {
    print(ii)
  }
  # print("iteration number")
  # print(ii)
  # cat("\n \n ")

}
# thinning  / warmup removing.

barplot(table(tail(unlist(d_accepted), n_accept - 500)))

plot(d_mat, type = "l")
plot(tail(unlist(d_accepted), n_accept - 5000), type = "l")


N <- 4
dang <- apply(returnSamplesOfLengthN(gamma_accepted_samples, N), 2, mean)

dang_full <- returnSamplesOfLengthN(gamma_accepted_samples, N)
dang_full_fits <- Q_full[, 1:N] %*% t(dang_full)

str(dang_full_fits)

dang_lower_fit <- apply(dang_full_fits, 1, function(x){quantile(x, 0.025)})
dang_upper_fit <- apply(dang_full_fits, 1, function(x){quantile(x, 0.975)})


dang_fit <- Q_full[, 1:N] %*% dang




dat_df$dang_fit <- dang_fit
dat_df$dang_lower_fit <- dang_lower_fit
dat_df$dang_upper_fit <- dang_upper_fit
ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = dang_fit))  +
  geom_ribbon(mapping = aes(x = x, ymin = dang_lower_fit, ymax = dang_upper_fit, alpha = 0.05))


# Traces
matplot(dang_full, type = "l")

#asdf
# plotting fits and ci's



mean_gamma <- apply(gamma_mat, 2, mean)
dat_df$y_fit_gamma <- Q %*% mean_gamma

ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = inital_fitteds))   +
  geom_line(mapping = aes(x = x, y = y_fit_gamma), colour = "red")


fitted_y_vals_gamma <- Q %*% t(gamma_mat)
dat_df$lower_int <- apply(fitted_y_vals_gamma, 1, function(x){quantile(x, 0.1)})
dat_df$upper_int <-apply(fitted_y_vals_gamma, 1, function(x){quantile(x, 0.9)})

ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = inital_fitteds))   +
  geom_line(mapping = aes(x = x, y = y_fit_gamma), colour = "red") +
  geom_ribbon(mapping = aes(x = x, ymin = lower_int, ymax = upper_int), alpha = 0.2)


# plotting traces, acfs, diagnostics

plot(gamma_mat[,d + 1], type = "l")
plot(var_mat[20000:n_iter], type = "l")


barplot(table(tail(unlist(d_accepted), n_accept - 200)))
