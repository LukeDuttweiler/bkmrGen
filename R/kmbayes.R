kmbayes <- function(y,
                    Z,
                    X,
                    K = 1,
                    iter = 1000,
                    warmup = 50,
                    nchains = 1,
                    family = 'gaussian',
                    id,
                    verbose = TRUE,
                    Znew = NULL,
                    starting.values = NULL,
                    control.params = NULL,
                    varsel = FALSE,
                    groups = NULL,
                    knots = NULL,
                    ztest = NULL,
                    rmethod = 'varying',
                    est.h = FALSE){
  argg <- as.list(match.call())
  ###################
  #Input Error Checks
  ###################

  ###################
  #Warnings
  ###################

  ###################
  #Split Sample
  ###################
  N <- length(y)

  groupSizes <- rep(floor(N/K), K)
  if(N %% K != 0){
    groupSizes[1:(N%%K)] <- groupSizes[1:(N%%K)] + 1
  }
  groupIdx <- cumsum(unlist(lapply(1:K, function(k){c(1, rep(0, groupSizes[k]-1))})))

  ###################
  #MCMC Sampling
  ###################
  #Use family name to specify sampler
  sampler <- eval(parse(text = paste0('bkmr_mcmc_', family)))

  samples <- lapply(1:K, function(k){
    #Prepare arguments for sampler (per group)
    sampCall <- argg[-1]
    sampCall$y <- y[groupIdx == k]
    sampCall$Z <- Z[groupIdx == k,]
    sampCall$X <- X[groupIdx == k,]
    sampCall$family <- NULL
    sampCall$K <- NULL

    return(do.call(sampler, sampCall))
  })

  #######################
  #Reconstruct Posterior
  #######################

  ###################
  #Format Return
  ###################
  return(samples)
}
