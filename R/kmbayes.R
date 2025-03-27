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
                    est.h = FALSE,
                    WASPSolver = 'lpsolve'){
  time1 <- Sys.time()
  argg <- as.list(match.call())
  ###################
  #Input Error Checks
  ###################

  #Check if ROI package exists, alert if not.
  if(K!= 1){
    ROIName <- paste0('ROI.plugin.', WASPSolver)
    if(!ROIName %in% rownames(installed.packages())){
      stop(paste0('In order to use optimizer ', WASPSolver,
                  ", please run install.packages('", ROIName,
                  "'), or pick a different solver."))
    }
  }
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

  #Get Beta Posteriors
  betaSamps <- lapply(samples, getElement, 'beta')

  betaPosts <- lapply(1:ncol(betaSamps[[1]]), function(i){
    sampI <- lapply(betaSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver)
    return(postI)
  })

  time2 <- Sys.time()

  ###################
  #Format Return
  ###################
  ret <- list('h.hat' = NULL,
              'beta' = NULL,
              'lambda' = NULL,
              'sigsq.eps' = NULL,
              'r' = NULL,
              'acc.r' = NULL,
              'acc.lambda' = NULL,
              'delta' = NULL,
              'acc.rdelta' = NULL,
              'move.type' = NULL,
              'est.h' = NULL,
              'time1' = time2-time1)
  ret <- c(ret, argg[-1])
  return(ret)
}
