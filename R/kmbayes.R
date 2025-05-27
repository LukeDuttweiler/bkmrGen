#' Fit Bayesian Kernel Machine Regression
#'
#' Main function to fit BKMR using MCMC methods and divide-and-conquer approaches with WASP
#'
#' TODO:
#' Deep dive on logit and poisson sampling functions
#' Update Traceplots and diagnostic stuff with genMCMCDiag
#' Implement multi-chain for original bkmr stuff
#' Implement warmup for original bkmr stuff
#' Make sure we can input different starting params or control values successfully (Look at the Validate___ functions)
#'
#' @param y a vector of outcome data of length \code{n}.
#' @param Z an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
#' @param X an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
#' @param K Number of splits to make in sample. Results are recombined using WASP.
#' @param n_samps Total number of samples to return from WASP. Only required if K > 1. Defaults to iter.
#' @param iter number of iterations to run the sampler
#' @param warmup Number of iterations to discard as warmups. Defaults to 100
#' @param nchains Number of MCMC or HMC chains to utilize. Defaults to 1.
#' @param family family object (see ?family for examples) with family and link information.
#' @param id optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
#' @param Znew optional matrix of new predictor values at which to predict \code{h}, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using \code{\link{SamplePred}}
#' @param starting.values list of starting values for each parameter. If not specified default values will be chosen.
#' @param control.params list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
#' @param varsel TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
#' @param groups optional vector (of length \code{M}) of group indicators for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed.
#' @param knots optional matrix of knot locations for implementing the Gaussian predictive process of Banerjee et al. (2008). Currently only implemented for models without a random intercept.
#' @param ztest optional vector indicating on which variables in Z to conduct variable selection (the remaining variables will be forced into the model).
#' @param rmethod for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
#' @param est.h TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
#' @param WASPSolver Linear-program solver to use. Solvers available are provided by the R package ROI (see ?ROI for details).
#'
#' @return an object of class "bkmrfit" (containing the posterior samples from the model fit), which has the associated methods: print, summary
#' @export
#'
#' @seealso For guided examples, go to \url{https://jenfb.github.io/bkmr/overview.html}
#' @references Bobb, JF, Valeri L, Claus Henn B, Christiani DC, Wright RO, Mazumdar M, Godleski JJ, Coull BA (2015). Bayesian Kernel Machine Regression for Estimating the Health Effects of Multi-Pollutant Mixtures. Biostatistics 16, no. 3: 493-508.
#' @references Banerjee S, Gelfand AE, Finley AO, Sang H (2008). Gaussian predictive process models for large spatial data sets. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 70(4), 825-848.
#'
#' @examples
#' ## First generate dataset
#' set.seed(111)
#' dat <- SimData(n = 50, M = 4)
#' y <- dat$y
#' Z <- dat$Z
#' X <- dat$X
#'
#'
#'
#' ## Fit model with component-wise variable selection
#' ## Using only 100 iterations to make example run quickly
#' ## Typically should use a large number of iterations for inference
#' set.seed(111)
#' fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 100, verbose = FALSE, varsel = TRUE)
kmbayes <- function(y,
                    Z,
                    X,
                    K = 1,
                    iter = 1000,
                    warmup = 100,
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
  argg$time1 <- Sys.time()


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
    stop("family ", family, " not yet implemented; must specify 'gaussian', 'binomial', or 'poisson'")
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

    if(link == 'logit' & varsel){#for a logit link, specify if using variable selection
      link <- 'logit_comp'
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

  ###################################
  #Prep Variables That May Be Needed
  ###################################
  argg$data.comps <- createDataComps(id = id, knots = knots)

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

  groupRuns <- lapply(1:K, function(k){
    #Prepare arguments for sampler (per group)
    sampCall <- argg
    sampCall$y <- y[groupIdx == k]
    sampCall$Z <- Z[groupIdx == k,, drop = F]
    sampCall$X <- X[groupIdx == k,, drop = F]
    sampCall$link <- link
    sampCall$missingX <- missingX

    return(do.call(sampler, sampCall))
  })

  #######################
  #Recover Full Control Params and Starting Values
  #######################
  argg$starting.values <- groupRuns[[1]]$starting.values
  argg$control.params <- groupRuns[[1]]$control.params

  #######################
  #Recover Samples
  #######################
  samples <- lapply(groupRuns, getElement, 'sampOutput')

  #Correct names to match STAN outputs to bkmr standard
  samples <- lapply(samples, function(kSet){
    names(kSet) <- gsub('_', '.', names(kSet))
    return(kSet)
  })

  #######################
  #Reconstruct Posterior
  #######################

  #Get Beta Posteriors
  betaSamps <- lapply(samples, getElement, 'beta')

  betaPosts <- sapply(1:ncol(betaSamps[[1]]), function(i){
    sampI <- lapply(betaSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver, n_samps = n_samps)
    return(postI)
  })

  #Get Lambda Posteriors
  lambdaSamps <- lapply(samples, getElement, 'lambda')
  lambdaSamps <- lapply(lambdaSamps, as.matrix)

  lambdaPosts <- sapply(1:ncol(lambdaSamps[[1]]), function(i){
    sampI <- lapply(lambdaSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver, n_samps = n_samps)
    return(postI)
  })

  #Get sigsq Posteriors
  sigsqSamps <- lapply(samples, getElement, 'sigsq.eps')
  sigsqSamps <- lapply(sigsqSamps, as.matrix)

  sigsqPosts <- sapply(1:ncol(sigsqSamps[[1]]), function(i){
    sampI <- lapply(sigsqSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver, n_samps = n_samps)
    return(postI)
  })

  #Get r Posteriors
  rSamps <- lapply(samples, getElement, 'r')
  rSamps <- lapply(rSamps, as.matrix)

  rPosts <- sapply(1:ncol(rSamps[[1]]), function(i){
    sampI <- lapply(rSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver, n_samps = n_samps)
    return(postI)
  })

  #Get delta Posteriors
  deltaSamps <- lapply(samples, getElement, 'delta')
  deltaSamps <- lapply(deltaSamps, as.matrix)

  deltaPosts <- sapply(1:ncol(deltaSamps[[1]]), function(i){
    sampI <- lapply(deltaSamps, function(b){return(b[,i])})
    postI <- wasp_univariate(sampI, solver = WASPSolver, n_samps = n_samps)
    return(postI)
  })

  #Get ystar Posteriors
  if(family != 'gaussian'){
    ystarSamps <- lapply(samples, getElement, 'ystar')
    ystarSamps <- lapply(ystarSamps, as.matrix)
    ystar <- do.call('cbind', ystarSamps)
  }else{
    ystar <- matrix(y, nrow = 1)
  }

  #Get h Posteriors
  #hSamps <- lapply(samples, getElement, 'h.hat')
  #hSamps <- lapply(hSamps, as.matrix)

  #Make sure samples is all matrices
  samples <- lapply(samples, function(subst){
    return(lapply(subst, as.matrix))
  })


  ###################
  #Format Return
  ###################
  ret <- list(#'h.hat' = hPosts,
              'ystar' = ystar,
              'beta' = betaPosts,
              'lambda' = lambdaPosts,
              'sigsq.eps' = sigsqPosts,
              'r' = rPosts,
              'delta' = deltaPosts,
              'time2' = Sys.time(),
              'subSamples' = samples)
  ret <- c(ret, argg)
  class(ret) <- c('bkmrfit', class(ret))
  return(ret)
}
