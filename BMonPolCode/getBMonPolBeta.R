library(rstan)

getBMonPolBeta <- function(x, y, q = d, iter) {
  N <- length(x)


  lower.bound <- 0
  upper.bound <- Inf

  alpha <- 1

  ##
  ## Given order of polynomial and bounds, figure out K and op.mode
  ## This is just following Proposition 2.1
  ##

  if (lower.bound == -Inf) {
    if (upper.bound == Inf) {
      op.mode <- 1 # whole real line
      k <- (q - 1)/2

    } else {
      op.mode <- 2 # (-inf, b]
      if(q %% 2 == 0) {
        k <- q / 2

      }  else {
        k <- (q - 1) / 2

      }
    }
  } else {
    if (upper.bound == Inf) {
      op.mode <- 3 # [a, Inf)
      if(q %% 2 == 0) {
        k <- q / 2

      }  else {
        k <- (q - 1) / 2

      }
    } else {
      op.mode <- 4 # [a,b]
      if(q %% 2 == 0) {
        k <- q / 2

      }  else {
        k <- (q - 1) / 2

      }
    }
  }

  if ( (q %% 2 == 0) & (op.mode == 1) ) {
    stop("Fitting a monotonic polynomial to the whole real line requires q to be odd")
  }

  ##
  ## Grid of values on which we want to have predictions
  ##

  Nnew <- 3
  xnew <- seq(from=min(x),to=max(x),length=Nnew)

  ##
  ## Put all the data to be passed to stan in a named list
  ##
  data.in <- list(N = N, q = q, K = k,
                  operation_mode = op.mode,
                  y = y, x = x,
                  a = lower.bound, b = upper.bound,
                  alpha = alpha, Nnew = Nnew,
                  xnew = xnew)

  ##
  ## Decide on number of iterations, this depends partly on the degree
  ## of the polynomial.
  ## See the 00README file
  ##



  ##
  ## Fit the model with stan
  ##
  ##
  ## Revert some of these global assigns onces debugging is completed
  prefit <<- stan("./BMonPolCode/MonPolyV0_0_4.stan", data = data.in, chains = 0)

  model.fit <- stan(fit = prefit, data = data.in,
                    iter = iter, seed = 170117, chains = 4, cores = 4,
                    save_warmup  = FALSE, warmup = 0.95*iter, refresh = 100,
                    control = list(max_treedepth = 15, adapt_delta = 0.9))


  ##
  ## beta final stuff, used for plotting

  beta.final.samples <- extract(model.fit, "beta_final")
  beta.final <- apply(beta.final.samples[[1]], 2, mean)

  return(beta.final)
}
