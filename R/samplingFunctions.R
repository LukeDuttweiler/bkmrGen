bkmr_mcmc_logit <- function(y,
                            Z,
                            X,
                            starting.values = list(rho = 2),
                            nchains = 1,
                            iter = 1000,
                            warmup = 50,
                            ...){

  ###################
  #Format for STAN
  ###################
  stanDat <- prepForStan(y, Z, X, starting.values)

  ###################
  #Run STAN
  ###################
  ft <- rstan::sampling(stanmodels$fit_logit, data = stanDat,
                        iter = iter + warmup, warmup = warmup,
                        chains = nchains, refresh = 0)

  ###################
  #Get MCMC Samples
  ###################
  samples <- rstan::extract(ft)
  return(samples)
}

bkmr_mcmc_poisson <- function(y,
                              Z,
                              X,
                              starting.values = list(rho = 2),
                              nchains = 1,
                              iter = 1000,
                              warmup = 50,
                              ...){
  ###################
  #Format for STAN
  ###################
  stanDat <- prepForStan(y, Z, X, starting.values)

  ###################
  #Run STAN
  ###################
  ft <- rstan::sampling(stanmodels$fit_poisson, data = stanDat,
                        iter = iter + warmup, warmup = warmup,
                        chains = nchains, refresh = 0)

  ###################
  #Get MCMC Samples
  ###################
  samples <- rstan::extract(ft)
  return(samples)
}
