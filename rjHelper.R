# functions specifically to enable the Reversible jump aspect of this routine

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

