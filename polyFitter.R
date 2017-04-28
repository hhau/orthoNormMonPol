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

source("basisGenerator.R")
source("mhHelpers.R")
source("plotDiagnoseFunctions.R")

# data thingos
x <- with(onechild, (day - min(day))/diff(range(day)))
y <- onechild$height
miny <- min(y)
ry <- diff(range(y))
y <- (y - miny)/ry
N <- length(x)

dat_df <- data.frame(x = x, y = y)
ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
  geom_line(mapping = aes(x = x, y = temp_fitted))



# mcmc params
n_iter <- 50000
d <- 8
innov_sigma <- 10 * .Machine$double.eps

#QR mats
qr_mats <- genAllMatrices(x, d)
Q <- qr_mats[[1]]
R_inv <- qr_mats[[2]]

gamma_mat <- matrix(NA, nrow = n_iter + 1, ncol = d + 1)
beta_mat <- matrix(NA, nrow = n_iter + 1, ncol = d + 1)

gamma_mat[1, ] <- genInitalGamma(dat_df, d, R_inv, Q)
beta_mat[1, ] <- gammaToBeta(gamma_mat[1, ], R_inv)

## mh loop
## there might need to be some dynamically changing innov stuff going on
## here
for (ii in 2:(n_iter + 1)) {
    gamma_proposal <- genMonoProp(gamma_mat[ii - 1, ], innov_sigma, R_inv)


}

