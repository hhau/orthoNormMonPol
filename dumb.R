# this is dumb, but i'm going to do it anyway
clearAndSetupSpace <- function() {
  # watch me global assign everything
  # Also clear everything before hand

  rm(list = ls())
  x <<- with(onechild, (day - min(day))/diff(range(day)))
  y <<- onechild$height
  miny <<- min(y)
  ry <<- diff(range(y))
  y <<- (y - miny)/ry
  N <<- length(x)

  dat_df <<- data.frame(x = x, y = y)
  ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y))
  #geom_line(mapping = aes(x = x, y = temp_fitted))


  # RJ params
  d_max <<- 12
  d_min <<- 1
  d <<- d_max


  # MH mcmc params
  n_iter <<- 500000

  # these need super careful selection?
  innov_sd_beta <<- 0.05
  innov_sd_var <<- 0.01

  ## inital value generation
  ## This takes fucking ages
  # beta_initial <- getBMonPolBeta(x = x, y = y, q = d, iter = 1000)

  # stan_rdump("beta_initial", file = "betainits.R")
  beta_initial <- read_rdump("./betainits.R")$beta_initial
  ## Faster but wrong?
  ## seems to want me to rescale my data to (-1, 1)
  # beta_initial_fast <<- monpol(y ~ x, data = as.data.frame(x = x ,y = y), degree = d_max,
                              # a = min(x), b = max(x))


  dat_df$initial_fitteds <<- evalPol(x, beta_initial)
  ggplot(data = dat_df) + geom_point(mapping = aes(x = x, y = y)) +
    geom_line(mapping = aes(x = x, y = initial_fitteds))  # +
  # geom_line(mapping = aes(x = x, y = Q %*% gamma_proposal), colour = "red")



  #QR mats
  qr_mats <<- genAllMatrices(x, d)
  Q_full <<- qr_mats[[1]]
  R_inv_full <<- qr_mats[[2]]

  ## inintal values
  # these are super important because we use them for independent proposals later
  gamma_inital <<- betaToGamma(beta_initial, R_inv_full)

  var_inital <<- 0.05


  #### RERUN FROM HERE TO BLANK DATA STORING OBJECTS

  # need to populate this `triangular` list


  # this is going to look really ugly.
  gamma_samples_list <<- list()
  gamma_samples_list[[1]] <<- gamma_inital


  gamma_accepted_samples <<- list()
  d_accepted <<- list()

  ## Don't think i can use these anymore (with varying dimension)
  # gamma_mat <<- matrix(NA, nrow = n_iter + 1, ncol = d + 1)
  # beta_mat <<- matrix(NA, nrow = n_iter + 1, ncol = d + 1)


  var_mat <<- matrix(NA, nrow = n_iter + 1, ncol = 1)
  d_mat <<- matrix(NA, nrow = n_iter + 1, ncol = 1)

  gamma_samples_list[[1]] <<- gamma_inital

  # beta_mat[1, ] <<- beta_initial
  var_mat[1] <<- var_inital
  d_mat[1] <<- d

  ## mh loop
  ## there might need to be some dynamically changing innov stuff going on
  ## here
  ##
  ## We are also choosing flat priors for our parameters (with the monotonicity
  ## caveat), and using a Random walk proposal for gamma(and by proxy beta)
  ## and an indepedent proposal for the variance

  n_accept <<- 0

  return(NULL)

}
