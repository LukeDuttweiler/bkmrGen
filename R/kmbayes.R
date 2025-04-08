#' Fit Bayesian Kernel Machine Regression
#'
#' Main function to fit BKMR using MCMC methods and divide-and-conquer approaches with WASP
#'
#' TODO: MUCH more work on function, much more work on documentation
#'
#' @param y Outcome
#' @param Z
#' @param X
#' @param K
#' @param n_samps
#' @param iter
#' @param warmup
#' @param nchains
#' @param family
#' @param id
#' @param verbose
#' @param Znew
#' @param starting.values
#' @param control.params
#' @param varsel
#' @param groups
#' @param knots
#' @param ztest
#' @param rmethod
#' @param est.h
#' @param WASPSolver
#'
#' @return
#' @export
#'
#' @examples
kmbayes <- function(y,
                    Z,
                    X,
                    K = 1,
                    iter = 1000,
                    warmup = 50,
                    nchains = 1,
                    n_samps = iter,
                    family = gaussian(),
                    id = NULL,
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

  #All arguments
  argg <- as.list(environment())


  ###################
  #Basic Defaults
  ###################

  missingX <- is.null(X) #Flag no X supplied
  if (missingX) X <- matrix(0, length(y), 1) #0s for covariates if not supplied
  hier_varsel <- !is.null(groups) #Hierarchical variable selection if groups supplied

  ###################
  #Input Error Checks
  ###################

  ##Argument checks: required arguments without defaults
  ##check vector/matrix sizes
  stopifnot (length(y) > 0, is.numeric(y), anyNA(y) == FALSE)
  if (!inherits(Z, "matrix"))  Z <- as.matrix(Z)
  stopifnot (is.numeric(Z), nrow(Z) == length(y), anyNA(Z) == FALSE)
  if (!inherits(X, "matrix"))  X <- as.matrix(X)
  stopifnot (is.numeric(X), nrow(X) == length(y), anyNA(X) == FALSE)


  ##Argument checks: Unique situations
  if (!is.null(id)) {
    stopifnot(length(id) == length(y), anyNA(id) == FALSE)
    if (!is.null(knots)) {
      message ("knots cannot be specified with id, resetting knots to null")
      knots<-NA
    }
  }

  if (!is.null(Znew)) {
    if (!inherits(Znew, "matrix"))  Znew <- as.matrix(Znew)
    stopifnot(is.numeric(Znew), ncol(Znew) == ncol(Z), anyNA(Znew) == FALSE)
  }

  if (!is.null(knots)) {
    if (!inherits(knots, "matrix"))  knots <- as.matrix(knots)
    stopifnot(is.numeric(knots), ncol(knots )== ncol(Z), anyNA(knots) == FALSE)
  }

  if (!is.null(groups)) {
    if (varsel == FALSE) {
      message ("groups should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(groups), length(groups) == ncol(Z), anyNA(groups) == FALSE)
    }
  }

  if (!is.null(ztest)) {
    if (varsel == FALSE) {
      message ("ztest should only be specified if varsel=TRUE, resetting varsel to TRUE")
      varsel <- TRUE
    } else {
      stopifnot(is.numeric(ztest), length(ztest) <= ncol(Z), anyNA(ztest) == FALSE, max(ztest) <= ncol(Z) )
    }
  }

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

  ##Argument check: for those with defaults, write message and reset to default if invalid
  if (iter < 1) {
    message ("invalid input for iter, resetting to default value 1000")
    nsamp <- 1000
  } else {
    nsamp <- iter
  }
  if (rmethod != "varying" & rmethod != "equal" & rmethod != "fixed") {
    message ("invalid value for rmethod, resetting to default varying")
    rmethod <- "varying"
  }
  if (verbose != FALSE & verbose != TRUE) {
    message ("invalid value for verbose, resetting to default FALSE")
    verbose <- FALSE
  }
  if (varsel != FALSE & varsel != TRUE) {
    message ("invalid value for varsel, resetting to default FALSE")
    varsel <- FALSE
  }

  ######################################
  #FAMILY CHECKS AND WARNINGS
  ######################################

  #Unpack options from family if given a function
  if('family' %in% class(family)){
    link <- family$link
    family <- family$family
  }

  if (!family %in% c("gaussian", "binomial", 'poisson')) {
    stop("family", family, "not yet implemented; must specify 'gaussian', 'binomial', or 'poisson'")
  }

  #Gaussian options
  if(family == 'gaussian'){
    if(!exists('link')){
      link <- 'identity'
    }

    if(!(link %in% c('identity'))){
      stop('Link ', link, ' not yet available for gaussian family. Must choose "identity" link.')
    }
  }

  #Binomial options
  if (family == "binomial") {
    if(!exists('link')){
      link <- 'logit'
      message('Logit link chosen for binomial family by default.')
    }

    if(!(link %in% c('logit', 'probit'))){
      stop('Link ', link, ' not yet available for binomial family. Must choose "logit" or "probit" link.')
    }

    if (!all(y %in% c(0, 1))) {
      stop("When family == 'binomial', y must be a vector containing only zeros and ones")
    }
  }

  #Poisson options
  if (family == "poisson") {
    if(!exists('link')){
      link <- 'log'
    }

    if(!(link %in% c('log'))){
      stop('Link ', link, ' not yet available for poisson family. Must choose "log" link.')
    }

    if (!is.integer(y) | any(y < 0)) {
      stop("When family == 'poisson', y must be a vector containing only non-negative integers")
    }
  }

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
    sampCall <- argg
    sampCall$y <- y[groupIdx == k]
    sampCall$Z <- Z[groupIdx == k,, drop = F]
    sampCall$X <- X[groupIdx == k,, drop = F]
    sampCall$family <- NULL
    sampCall$K <- NULL
    sampCall$missingX <- missingX

    return(do.call(sampler, sampCall))
  })

  #######################
  #Reconstruct Posterior
  #######################

  #Get Beta Posteriors
  betaSamps <- lapply(samples, getElement, 'beta')

  betaPosts <- lapply(1:ncol(betaSamps[[1]]), function(i){
    sampI <- lapply(betaSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver, n_samps = n_samps)
    return(postI)
  })

  time2 <- Sys.time()

  ###################
  #Format Return
  ###################
  ret <- list('h.hat' = NULL,
              'beta' = betaPosts,
              'lambda' = NULL,
              'sigsq.eps' = NULL,
              'r' = NULL,
              'acc.r' = NULL,
              'acc.lambda' = NULL,
              'delta' = NULL,
              'acc.rdelta' = NULL,
              'move.type' = NULL,
              'est.h' = NULL,
              'allSamples' = samples)
  ret <- c(ret, argg)
  return(ret)

  #########################################################
  #ORIGINAL PACKAGE CODE
  #########################################################


  control.params$r.params <- NULL
  chain$time2 <- Sys.time()
  chain$iter <- nsamp
  chain$family <- family
  chain$starting.values <- starting.values
  chain$control.params <- control.params
  chain$X <- X
  chain$Z <- Z
  chain$y <- y
  chain$ztest <- ztest
  chain$data.comps <- data.comps
  if (!is.null(Znew)) chain$Znew <- Znew
  if (!is.null(groups)) chain$groups <- groups
  chain$varsel <- varsel
  class(chain) <- c("bkmrfit", class(chain))
  chain
}
