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

source("basisGenerator.R")
source("mhHelpers.R")
source("plotDiagnoseFunctions.R")
source("rjHelper.R")

source("./BMonPolCode/getBMonPolBeta.R")


# data thingos
x <- with(onechild, (day - min(day))/diff(range(day)))
y <- onechild$height
miny <- min(y)
ry <- diff(range(y))
y <- (y - miny)/ry
N <- length(x)

dat_df <- data.frame(x = x, y = y)
ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y))
  #geom_line(mapping = aes(x = x, y = temp_fitted))


# RJ params
d_max <- 12
d_min <- 1
d <- d_max


# MH mcmc params
n_iter <- 100000

# these need super careful selection?
innov_sd_beta <- 0.05
innov_sd_var <- 0.01

## inital value generation
## This takes fucking ages
beta_initial <- getBMonPolBeta(x = x, y = y, q = d, iter = 1000)

stan_rdump("beta_initial", file = "betainits.R")
beta_initial <- read_rdump("./betainits.R")$beta_initial
## Faster but wrong?
## seems to want me to rescale my data to (-1, 1)
# beta_initial_fast <- monpol(y ~ x, data = as.data.frame(x = x ,y = y), degree = d_max,
#                            a = min(x), b = max(x))$beta


dat_df$initial_fitteds <- evalPol(x, beta_initial)
ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = initial_fitteds))  # +
  # geom_line(mapping = aes(x = x, y = Q %*% gamma_proposal), colour = "red")



#QR mats
qr_mats <- genAllMatrices(x, d)
Q_full <- qr_mats[[1]]
R_inv_full <- qr_mats[[2]]

gamma_inital <- betaToGamma(beta_initial, R_inv_full)

var_inital <- 0.05
# need to choose some gamma proposal values.
#

# need to populate this `triangular` list
most_recent_gamma_list <- list()

for (qq in 1:(d_max - d_min + 1)) {
  most_recent_gamma_list[[qq]] <- gamma_inital[1:(d_min + qq)]
}

# this is going to look really ugly.
gamma_samples_list <- list()
gamma_samples_list[[1]] <- gamma_inital

gamma_accepted_samples <- list()
d_accepted <- list()

## Don't think i can use these anymore (with varying dimension)
# gamma_mat <- matrix(NA, nrow = n_iter + 1, ncol = d + 1)
# beta_mat <- matrix(NA, nrow = n_iter + 1, ncol = d + 1)


var_mat <- matrix(NA, nrow = n_iter + 1, ncol = 1)
d_mat <- matrix(NA, nrow = n_iter + 1, ncol = 1)

gamma_samples_list[[1]] <- gamma_inital
most_recent_gamma_list[[(d - d_min + 1)]] <- gamma_inital

# beta_mat[1, ] <- beta_initial
var_mat[1] <- var_inital
d_mat[1] <- d

## mh loop
## there might need to be some dynamically changing innov stuff going on
## here
##
## We are also choosing flat priors for our parameters (with the monotonicity
## caveat), and using a Random walk proposal for gamma(and by proxy beta)
## and an indepedent proposal for the variance

n_accept <- 0
for (ii in 2:(n_iter + 1)) {

    d <- d_mat[ii - 1]
    d_prop <- dimProposer(d, d_min, d_max)

    variance_proposal <- genVarProp(var_prev = var_mat[ii - 1],
                                    innov_sd_var = innov_sd_var)

    if (d_prop < d) {
      # If i drop a dimension, take the mean as d_prop of the gammas before
      # for the moment I am assuming that q_{d_{k -> k'}} (u) is 1 when
      # k' is smaller than k, because u isn't zero, it doesn't exist
      #
      # Now to fix this horrific assumption.
      #
      # THIS IS WRONG RIGHT NOW, NEEDS TO USE SAMPLE LIST
      # ALSO VERY INCORRECT IN TERMS OF HOW THIS IS ALL CALCULATED.
      proposal_mean <- gamma_samples_list[[ii - 1]][1:(d_prop + 1)]

      gamma_proposal <- genMonoProp(gamma_old = proposal_mean,
                                    innov_sigma = innov_sd_beta,
                                    R_inv = R_inv_full[1:(d_prop + 1), 1:(d_prop + 1)])

      # now switch minds, the proposal is no longer the proposal, it is the old
      # value. We can invert the last step.

      # q_dim_swap now pertains to the probability of the last parameter being
      # generated?????????? the last parameter no longer exists

      # I need to find the sample before ii - 1 when we were last in dimension
      # d

      zz <- ii - 2
      if (ii > 50) {
        while (zz > 0) {
          if (length(gamma_samples_list[[zz]]) == (d + 1)) {
            that_fucking_mean_vec <- gamma_samples_list[[zz]]
            break
          } else {
            print("shitfuck")
            that_fucking_mean_vec <- gamma_samples_list[[ii - 1]]
          }
          zz <- zz - 1
        }
      } else {
        that_fucking_mean_vec <- gamma_samples_list[[ii - 1]]
      }

      q_dim_swap <- dnorm(x = gamma_samples_list[[ii - 1]][d + 1],
                          mean = that_fucking_mean_vec[d + 1],
                          sd = innov_sd_beta)

      accept_propability <- (calcAcceptProb(Q_old = Q_full[,1:(d_prop + 1)],
                                           Q_prop = Q_full[, 1:(d + 1)],
                                           gamma_prop = gamma_samples_list[[ii - 1]],
                                           gamma_old = gamma_proposal,
                                           var_prop = var_mat[ii - 1],
                                           var_old = variance_proposal,
                                           y_vec = y,
                                           q_dim_swap = q_dim_swap)) ^ (-1)

      accept_propability <- min(accept_propability,  1)

    } else if (d_prop > d) {
      # If i increase a dimension, the mean is then the same as where I was
      # before for the first d+1 components, and
      proposal_mean_base <- gamma_samples_list[[ii - 1]]
      proposal_mean_extra <- most_recent_gamma_list[[d_prop - d_min + 1]][d_prop + 1]

      proposal_mean <- c(proposal_mean_base, proposal_mean_extra)

      gamma_proposal <- genMonoProp(gamma_old = proposal_mean, innov_sigma = innov_sd_beta,
                                    R_inv = R_inv_full[1:(d_prop + 1), 1:(d_prop + 1)])

      q_dim_swap <- dnorm(x = gamma_proposal[d_prop + 1], mean = proposal_mean_extra,
                          sd = innov_sd_beta)

      accept_propability <- calcAcceptProb(Q_old = Q_full[,1:(d + 1)],
                                           Q_prop = Q_full[, 1:(d_prop + 1)],
                                           gamma_prop = gamma_proposal,
                                           gamma_old = gamma_samples_list[[ii - 1]],
                                           var_prop = variance_proposal,
                                           var_old = var_mat[ii - 1],
                                           y_vec = y,
                                           q_dim_swap = q_dim_swap)
      accept_propability <- min(accept_propability,  1)

    } else {
      # No dimension switch
      # business as usual
      proposal_mean <- gamma_samples_list[[ii - 1]]
      gamma_proposal <- genMonoProp(gamma_old = proposal_mean,
                                    innov_sigma = innov_sd_beta,
                                    R_inv = R_inv_full[1:(d_prop + 1), 1:(d_prop + 1)])

      q_dim_swap <- 1

      accept_propability <- calcAcceptProb(Q_old = Q_full[,1:(d + 1)],
                                           Q_prop = Q_full[, 1:(d_prop + 1)],
                                           gamma_prop = gamma_proposal,
                                           gamma_old = gamma_samples_list[[ii - 1]],
                                           var_prop = variance_proposal,
                                           var_old = var_mat[ii - 1],
                                           y_vec = y,
                                           q_dim_swap = q_dim_swap)
      accept_propability <- min(accept_propability,  1)
    }

    # somewhere in here I'll need to write a jump_probability step.
    # gamma_proposal <- genMonoProp(gamma_old = gamma_mat[ii - 1, ],
    #                               innov_sigma = innov_sd_beta,
    #                               R_inv = R_inv)

    #beta_propsal <- gammaToBeta(gamma_proposal, R_inv)









  if (runif(1) < accept_propability) {
    # gamma_mat[ii, ] <- gamma_proposal

    gamma_samples_list[[ii]] <- gamma_proposal
    most_recent_gamma_list[[d_prop - d_min + 1]] <- gamma_proposal



    d_mat[ii] <- d_prop


    var_mat[ii] <- variance_proposal
    n_accept <- n_accept + 1

    gamma_accepted_samples[[n_accept]] <- gamma_proposal
    d_accepted[[n_accept]] <- d_prop
  } else {
    # gamma_mat[ii, ] <- gamma_mat[ii - 1, ]
    var_mat[ii] <- var_mat[ii - 1]

    gamma_samples_list[[ii]] <- gamma_samples_list[[ii - 1]]

    d_mat[ii] <- d_mat[ii - 1]

  }

  # if (ii %% 1000 == 0) {
    print(ii)
  # }
}
# thinning  / warmup removing.

barplot(table(tail(unlist(d_accepted), n_accept - 5000)))

# plot(d_mat, type = "l")
plot(tail(unlist(d_accepted), n_accept - 5000), type = "l")


N <- 10
dang <- apply(returnSamplesOfLengthN(gamma_accepted_samples, N), 2, mean)

dang_fit <- Q_full[, 1:N] %*% dang

dat_df$dang_fit <- dang_fit
ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = dang_fit))  # +

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
