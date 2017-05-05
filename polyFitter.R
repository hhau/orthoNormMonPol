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



# MH mcmc params
n_iter <- 100000
d <- 8

# these need super careful selection?
innov_sd_beta <- 0.05
innov_sd_var <- 0.01

## inital value generation
beta_initial <- getBMonPolBeta(x = x, y = y, q = d, iter = 1000)
dat_df$inital_fitteds <- evalPol(x, beta_initial)
ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = inital_fitteds))  # +
  # geom_line(mapping = aes(x = x, y = Q %*% gamma_proposal), colour = "red")



#QR mats
qr_mats <- genAllMatrices(x, d)
Q <- qr_mats[[1]]
R_inv <- qr_mats[[2]]

gamma_inital <- betaToGamma(beta_initial, R_inv)

var_inital <- 0.05
# need to choose some gamma proposal values.
#


gamma_mat <- matrix(NA, nrow = n_iter + 1, ncol = d + 1)
beta_mat <- matrix(NA, nrow = n_iter + 1, ncol = d + 1)
var_mat <- matrix(NA, nrow = n_iter + 1, ncol = 1)

gamma_mat[1, ] <- gamma_inital
beta_mat[1, ] <- beta_initial
var_mat[1] <- var_inital

## mh loop
## there might need to be some dynamically changing innov stuff going on
## here
##
## We are also choosing flat priors for our parameters (with the monotonicity
## caveat), and using a Random walk proposal for gamma(and by proxy beta)
## and an indepedent proposal for the variance

n_accept <- 0
for (ii in 2:(n_iter + 1)) {

    # somewhere in here I'll need to write a jump_probability step.
    gamma_proposal <- genMonoProp(gamma_old = gamma_mat[ii - 1, ],
                                  innov_sigma = innov_sd_beta,
                                  R_inv = R_inv)

    #beta_propsal <- gammaToBeta(gamma_proposal, R_inv)

    variance_proposal <- genVarProp(var_prev = var_mat[ii - 1],
                                    innov_sd_var = innov_sd_var)

    accept_propability <- calcLikelihoodRatio(Q = Q, gamma_prop = gamma_proposal,
                                              gamma_old = gamma_mat[ii - 1,],
                                              var_prop = variance_proposal,
                                              var_old = var_mat[ii - 1],
                                              y_vec = y)
  if (runif(1) < accept_propability) {
    gamma_mat[ii, ] <- gamma_proposal
    var_mat[ii] <- variance_proposal
    n_accept <- n_accept + 1
  } else {
    gamma_mat[ii, ] <- gamma_mat[ii - 1, ]
    var_mat[ii] <- var_mat[ii - 1]
  }

  if (ii %% 1000 == 0) {
    print(ii)
  }
}
# thinning  / warmup removing.


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
