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
