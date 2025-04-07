bkmr_mcmc_gaussian <- function(){

}

bkmr_mcmc_probit <- function(){

}

#' BKMR Sampling on Logit Model using Hamiltonian Monte Carlo
#'
#' Fits bkmr_logit model using STAN. See file inst/stan/fit_logit.stan for STAN specifications.
#'
#' @inheritParams kmbayes
#'
#' @param ... Catch all to allow package flexibility
#'
#' @return STAN output from bkmr_logit fit
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

#' BKMR Sampling on Poisson Model using Hamiltonian Monte Carlo
#'
#' Fits bkmr_poisson model using STAN. See file inst/stan/fit_poisson.stan for STAN specifications.
#'
#' @inheritParams kmbayes
#'
#' @param ... Catch all to allow package flexibility
#'
#' @return STAN output from bkmr_poisson fit
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
