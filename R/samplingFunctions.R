#' MCMC code to fit Gaussian BKMR Model
#'
#' Code slightly modified from original bkmr package, Bobb 2015
#'
#' @param missingX Logical Flag to indicate if matrix X was provided to the function
#' @param data.comps Components needed for V matrix
#' @param ... Catches extra arguments
#' @inheritParams kmbayes
#'
#' @return Samples from gaussian MCMC bkmr sampler.
bkmr_mcmc_gaussian <- function(y,
                               Z,
                               X,
                               starting.values,
                               nchains,
                               iter,
                               warmup,
                               id,
                               varsel,
                               Znew,
                               groups,
                               rmethod,
                               est.h,
                               knots,
                               ztest,
                               control.params,
                               verbose,
                               missingX,
                               data.comps,
                               ...){

  ## start JB code

  ## create empty matrices to store the posterior draws in
  chain <- list(h.hat = matrix(0, iter, nrow(Z)),
                beta = matrix(0, iter, ncol(X)),
                lambda = matrix(NA, iter, data.comps$nlambda),
                sigsq.eps = rep(NA, iter),
                r = matrix(NA, iter, ncol(Z)),
                acc.r = matrix(0, iter, ncol(Z)),
                acc.lambda = matrix(0, iter, data.comps$nlambda),
                delta = matrix(1, iter, ncol(Z))
  )
  if (varsel) {
    chain$acc.rdelta <- rep(0, iter)
    chain$move.type <- rep(0, iter)
  }

  ## components to predict h(Znew)
  if (!is.null(Znew)) {
    if (inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
    chain$hnew <- matrix(0,iter,nrow(Znew))
    colnames(chain$hnew) <- rownames(Znew)
  }

  ## components if model selection is being done
  if (varsel) {
    if (is.null(ztest)) {
      ztest <- 1:ncol(Z)
    }
    rdelta.update <- rdelta.comp.update
  } else {
    ztest <- NULL
  }

  ## control parameters (lambda.jump default lower for probit model to improve mixing)
  control.params.default <- list(lambda.jump = rep(10, data.comps$nlambda),
                                 mu.lambda = rep(10, data.comps$nlambda),
                                 sigma.lambda = rep(10, data.comps$nlambda),
                                 a.p0 = 1, b.p0 = 1, r.prior = "invunif", a.sigsq = 1e-3,
                                 b.sigsq = 1e-3, mu.r = 5, sigma.r = 5, r.muprop = 1,
                                 r.jump = 0.1, r.jump1 = 2, r.jump2 = 0.1, r.a = 0, r.b = 100)
  if (!is.null(control.params)){
    control.params <- modifyList(control.params.default, as.list(control.params))
    validateControlParams(varsel, family, id, control.params)
  } else {
    control.params <- control.params.default
  }

  control.params$r.params <- with(control.params, list(mu.r = mu.r, sigma.r = sigma.r, r.muprop = r.muprop, r.jump = r.jump, r.jump1 = r.jump1, r.jump2 = r.jump2, r.a = r.a, r.b = r.b))

  ## components if grouped model selection is being done
  if (!is.null(groups)) {
    hier_varsel <- TRUE
    rdelta.update <- rdelta.group.update
    control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
  }else{
    hier_varsel <- FALSE
  }

  ## specify functions for doing the Metropolis-Hastings steps to update r
  e <- environment()
  rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
  rprior.logdens <- rfn$rprior.logdens
  environment(rprior.logdens) <- e
  rprop.gen1 <- rfn$rprop.gen1
  environment(rprop.gen1) <- e
  rprop.logdens1 <- rfn$rprop.logdens1
  environment(rprop.logdens1) <- e
  rprop.gen2 <- rfn$rprop.gen2
  environment(rprop.gen2) <- e
  rprop.logdens2 <- rfn$rprop.logdens2
  environment(rprop.logdens2) <- e
  rprop.gen <- rfn$rprop.gen
  environment(rprop.gen) <- e
  rprop.logdens <- rfn$rprop.logdens
  environment(rprop.logdens) <- e
  rm(e, rfn)

  ## initial values
  starting.values0 <- list(h.hat = 1, beta = NULL, sigsq.eps = NULL, r = 1, lambda = 10, delta = 1)
  if (is.null(starting.values)) {
    starting.values <- starting.values0
  } else {
    starting.values <- modifyList(starting.values0, starting.values)
    validateStartingValues(varsel, y, X, Z, starting.values, rmethod)
  }

  if (is.null(starting.values$beta) | is.null(starting.values$sigsq.eps)) {
    lmfit0 <- lm(y ~ Z + X)
    if (is.null(starting.values$beta)) {
      coefX <- coef(lmfit0)[grep("X", names(coef(lmfit0)))]
      starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
    }
    if (is.null(starting.values$sigsq.eps)) {
      starting.values$sigsq.eps <- summary(lmfit0)$sigma^2
    }
  }

  ##truncate vectors that are too long
  if (length(starting.values$h.hat) > length(y)) {
    starting.values$h.hat <- starting.values$h.hat[1:length(y)]
  }
  if (length(starting.values$beta) > ncol(X)) {
    starting.values$beta <- starting.values$beta[1:ncol(X)]
  }
  if (length(starting.values$delta) > ncol(Z)) {
    starting.values$delta <- starting.values$delta[1:ncol(Z)]
  }
  if (varsel==FALSE & rmethod == "equal" & length(starting.values$r) > 1) {
    starting.values$r <- starting.values$r[1] ## this should only happen if rmethod == "equal"
  } else if (length(starting.values$r) > ncol(Z)) {
    starting.values$r <- starting.values$r[1:ncol(Z)]
  }

  chain$h.hat[1, ] <- starting.values$h.hat
  chain$beta[1, ] <- starting.values$beta
  chain$lambda[1, ] <- starting.values$lambda
  chain$sigsq.eps[1] <- starting.values$sigsq.eps
  chain$r[1, ] <- starting.values$r
  if (varsel) {
    chain$delta[1,ztest] <- starting.values$delta
  }
  if (!is.null(groups)) {
    ## make sure starting values are consistent with structure of model
    if (!all(sapply(unique(groups), function(x) sum(chain$delta[1, ztest][groups == x])) == 1)) {
      # warning("Specified starting values for delta not consistent with model; using default")
      starting.values$delta <- rep(0, length(groups))
      starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
    }
    chain$delta[1,ztest] <- starting.values$delta
    chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
  }
  chain$est.h <- est.h

  ## components
  Vcomps <- makeVcomps(r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps)

  # set print progress options
  opts <- set_verbose_opts(
    verbose_freq = control.params$verbose_freq,
    verbose_digits = control.params$verbose_digits,
    verbose_show_ests = control.params$verbose_show_ests,
    tot_iter=iter
  )

  ## start sampling ####
  for (s in 2:iter) {

    ## continuous version of outcome
    ycont <- y

    ## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)

    ## beta
    if (!missingX) {
      chain$beta[s,] <- beta.update(X = X, Vinv = Vcomps$Vinv, y = ycont, sigsq.eps = chain$sigsq.eps[s - 1])
    }

    ## \sigma_\epsilon^2
    chain$sigsq.eps[s] <- sigsq.eps.update(y = ycont, X = X, beta = chain$beta[s,], Vinv = Vcomps$Vinv, a.eps = control.params$a.sigsq, b.eps = control.params$b.sigsq)


    ## lambda
    lambdaSim <- chain$lambda[s - 1,]
    for (comp in 1:data.comps$nlambda) {
      varcomps <- lambda.update(r = chain$r[s - 1,], delta = chain$delta[s - 1,], lambda = lambdaSim, whichcomp = comp, y = ycont, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params)
      lambdaSim <- varcomps$lambda
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.lambda[s,comp] <- varcomps$acc
      }
    }
    chain$lambda[s,] <- lambdaSim

    ## r
    rSim <- chain$r[s - 1,]
    comp <- which(!1:ncol(Z) %in% ztest)
    if (length(comp) != 0) {
      if (rmethod == "equal") { ## common r for those variables not being selected
        varcomps <- r.update(r = rSim, whichcomp = comp, delta = chain$delta[s - 1,],
                             lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,],
                             sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z,
                             data.comps = data.comps, control.params = control.params,
                             rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
                             rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
                             rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
                             rprop.logdens = rprop.logdens)
        rSim <- varcomps$r
        if (varcomps$acc) {
          Vcomps <- varcomps$Vcomps
          chain$acc.r[s, comp] <- varcomps$acc
        }
      } else if (rmethod == "varying") { ## allow a different r_m
        for (whichr in comp) {
          varcomps <- r.update(r = rSim, whichcomp = whichr, delta = chain$delta[s - 1,],
                               lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,],
                               sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z,
                               data.comps = data.comps, control.params = control.params,
                               rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
                               rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
                               rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
                               rprop.logdens = rprop.logdens)
          rSim <- varcomps$r
          if (varcomps$acc) {
            Vcomps <- varcomps$Vcomps
            chain$acc.r[s, whichr] <- varcomps$acc
          }
        }
      }
    }
    ## for those variables being selected: joint posterior of (r,delta)
    if (varsel) {
      varcomps <- rdelta.update(r = rSim, delta = chain$delta[s - 1,],
                                lambda = chain$lambda[s,], y = ycont, X = X,
                                beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s],
                                Vcomps = Vcomps, Z = Z, ztest = ztest,
                                data.comps = data.comps, control.params = control.params,
                                rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
                                rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
                                rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
                                rprop.logdens = rprop.logdens)
      chain$delta[s,] <- varcomps$delta
      rSim <- varcomps$r
      chain$move.type[s] <- varcomps$move.type
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.rdelta[s] <- varcomps$acc
      }
    }
    chain$r[s,] <- rSim

    ###################################################
    ## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)

    if (est.h) {
      hcomps <- h.update(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s], y = ycont, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z, data.comps = data.comps)
      chain$h.hat[s,] <- hcomps$hsamp
      if (!is.null(hcomps$hsamp.star)) { ## GPP
        Vcomps$hsamp.star <- hcomps$hsamp.star
      }
      rm(hcomps)
    }

    ###################################################
    ## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)

    if (!is.null(Znew)) {
      chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s], r = chain$r[s,], y = ycont, X = X, beta = chain$beta[s,], data.comps = data.comps)
    }

    ###################################################
    ## print details of the model fit so far
    print_diagnostics(verbose = verbose, opts = opts, curr_iter = s, tot_iter = iter, chain = chain, varsel = varsel, hier_varsel = hier_varsel, ztest = ztest, Z = Z, groups = groups)

  }

  #Remove r.params from control.params for later functions
  control.params$r.params <- NULL

  #Put data into chain return (so subsets are easily identified)
  chain$y <- y
  chain$X <- X
  chain$Z <- Z

  return(list(sampOutput = chain,
              starting.values = starting.values,
              control.params = control.params))
}

#' Redirect Binomial family using specified link function (default will be 'logit')
#'
#' @param link character vector giving link function name. Currently implemented for 'logit' and 'probit'
#' @param ... All arguments passed through to sampler function
#'
#' @return Output from sampler function
bkmr_mcmc_binomial <- function(link, ...){
  #Redirect based on link function
  samplerLink <- eval(parse(text = paste0('bkmr_mcmc_', link)))

  return(samplerLink(...))
}

#' MCMC code to fit Probit BKMR Model
#'
#' Code slightly modified from original bkmr package, Bobb 2015
#'
#' @param missingX Logical Flag to indicate if matrix X was provided to the function
#' @param data.comps Components needed for V matrix
#' @param ... Catches extra arguments
#' @inheritParams kmbayes
#'
#' @return Samples from MCMC probit bkmr model
bkmr_mcmc_probit <- function(y,
                             Z,
                             X,
                             starting.values,
                             nchains,
                             iter,
                             warmup,
                             id,
                             varsel,
                             Znew,
                             groups,
                             rmethod,
                             est.h,
                             knots,
                             ztest,
                             control.params,
                             verbose,
                             missingX,
                             data.comps,
                             ...){

  ## start JB code

  ## create empty matrices to store the posterior draws in
  chain <- list(h.hat = matrix(0, iter, nrow(Z)),
                beta = matrix(0, iter, ncol(X)),
                lambda = matrix(NA, iter, data.comps$nlambda),
                sigsq.eps = rep(NA, iter),
                r = matrix(NA, iter, ncol(Z)),
                acc.r = matrix(0, iter, ncol(Z)),
                acc.lambda = matrix(0, iter, data.comps$nlambda),
                delta = matrix(1, iter, ncol(Z))
  )
  if (varsel) {
    chain$acc.rdelta <- rep(0, iter)
    chain$move.type <- rep(0, iter)
  }

  #Set up latent variable since this is probit model
  chain$ystar <- matrix(0, iter, length(y))

  ## components to predict h(Znew)
  if (!is.null(Znew)) {
    if (inherits(Znew, "data.frame")) Znew <- data.matrix(Znew)
    chain$hnew <- matrix(0,iter,nrow(Znew))
    colnames(chain$hnew) <- rownames(Znew)
  }

  ## components if model selection is being done
  if (varsel) {
    if (is.null(ztest)) {
      ztest <- 1:ncol(Z)
    }
    rdelta.update <- rdelta.comp.update
  } else {
    ztest <- NULL
  }

  ## control parameters (lambda.jump default lower for probit model to improve mixing)
  control.params.default <- list(lambda.jump = rep(sqrt(10), data.comps$nlambda),
                                 mu.lambda = rep(10, data.comps$nlambda),
                                 sigma.lambda = rep(10, data.comps$nlambda),
                                 a.p0 = 1, b.p0 = 1, r.prior = "invunif", a.sigsq = 1e-3,
                                 b.sigsq = 1e-3, mu.r = 5, sigma.r = 5, r.muprop = 1,
                                 r.jump = 0.1, r.jump1 = 2, r.jump2 = 0.1, r.a = 0, r.b = 100)
  if (!is.null(control.params)){
    control.params <- modifyList(control.params.default, as.list(control.params))
    validateControlParams(varsel, family, id, control.params)
  } else {
    control.params <- control.params.default
  }

  control.params$r.params <- with(control.params, list(mu.r = mu.r, sigma.r = sigma.r, r.muprop = r.muprop, r.jump = r.jump, r.jump1 = r.jump1, r.jump2 = r.jump2, r.a = r.a, r.b = r.b))

  ## components if grouped model selection is being done
  if (!is.null(groups)) {
    hier_varsel <- TRUE
    rdelta.update <- rdelta.group.update
    control.params$group.params <- list(groups = groups, sel.groups = sapply(unique(groups), function(x) min(seq_along(groups)[groups == x])), neach.group = sapply(unique(groups), function(x) sum(groups %in% x)))
  }else{
    hier_varsel <- FALSE
  }

  ## specify functions for doing the Metropolis-Hastings steps to update r
  e <- environment()
  rfn <- set.r.MH.functions(r.prior = control.params$r.prior)
  rprior.logdens <- rfn$rprior.logdens
  environment(rprior.logdens) <- e
  rprop.gen1 <- rfn$rprop.gen1
  environment(rprop.gen1) <- e
  rprop.logdens1 <- rfn$rprop.logdens1
  environment(rprop.logdens1) <- e
  rprop.gen2 <- rfn$rprop.gen2
  environment(rprop.gen2) <- e
  rprop.logdens2 <- rfn$rprop.logdens2
  environment(rprop.logdens2) <- e
  rprop.gen <- rfn$rprop.gen
  environment(rprop.gen) <- e
  rprop.logdens <- rfn$rprop.logdens
  environment(rprop.logdens) <- e
  rm(e, rfn)

  ## initial values
  starting.values0 <- list(h.hat = 1, beta = NULL, sigsq.eps = NULL, r = 1, lambda = 10, delta = 1)
  if (is.null(starting.values)) {
    starting.values <- starting.values0
  } else {
    starting.values <- modifyList(starting.values0, starting.values)
    validateStartingValues(varsel, y, X, Z, starting.values, rmethod)
  }

  #Starting values for probit model
  starting.values$sigsq.eps <- 1 ## always equal to 1
  if (is.null(starting.values$beta) | is.null(starting.values$ystar)) {
    probitfit0 <- try(glm(y ~ Z + X, family = binomial(link = "probit")))
    if (!inherits(probitfit0, "try-error")) {
      if (is.null(starting.values$beta)) {
        coefX <- coef(probitfit0)[grep("X", names(coef(probitfit0)))]
        starting.values$beta <- unname(ifelse(is.na(coefX), 0, coefX))
      }
      if (is.null(starting.values$ystar)) {
        starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
      }
    } else {
      starting.values$beta <- 0
      starting.values$ystar <- ifelse(y == 1, 1/2, -1/2)
    }
  }

  ##truncate vectors that are too long
  if (length(starting.values$h.hat) > length(y)) {
    starting.values$h.hat <- starting.values$h.hat[1:length(y)]
  }
  if (length(starting.values$beta) > ncol(X)) {
    starting.values$beta <- starting.values$beta[1:ncol(X)]
  }
  if (length(starting.values$delta) > ncol(Z)) {
    starting.values$delta <- starting.values$delta[1:ncol(Z)]
  }
  if (varsel==FALSE & rmethod == "equal" & length(starting.values$r) > 1) {
    starting.values$r <- starting.values$r[1] ## this should only happen if rmethod == "equal"
  } else if (length(starting.values$r) > ncol(Z)) {
    starting.values$r <- starting.values$r[1:ncol(Z)]
  }

  chain$h.hat[1, ] <- starting.values$h.hat
  chain$beta[1, ] <- starting.values$beta
  chain$lambda[1, ] <- starting.values$lambda
  chain$sigsq.eps[1] <- starting.values$sigsq.eps
  chain$r[1, ] <- starting.values$r
  if (varsel) {
    chain$delta[1,ztest] <- starting.values$delta
  }

  #Starting values for probit model specifically
  chain$ystar[1, ] <- starting.values$ystar
  chain$sigsq.eps[] <- starting.values$sigsq.eps ## does not get updated

  if (!is.null(groups)) {
    ## make sure starting values are consistent with structure of model
    if (!all(sapply(unique(groups), function(x) sum(chain$delta[1, ztest][groups == x])) == 1)) {
      # warning("Specified starting values for delta not consistent with model; using default")
      starting.values$delta <- rep(0, length(groups))
      starting.values$delta[sapply(unique(groups), function(x) min(which(groups == x)))] <- 1
    }
    chain$delta[1,ztest] <- starting.values$delta
    chain$r[1,ztest] <- ifelse(chain$delta[1,ztest] == 1, chain$r[1,ztest], 0)
  }
  chain$est.h <- est.h

  ## components
  Vcomps <- makeVcomps(r = chain$r[1, ], lambda = chain$lambda[1, ], Z = Z, data.comps = data.comps)

  # set print progress options
  opts <- set_verbose_opts(
    verbose_freq = control.params$verbose_freq,
    verbose_digits = control.params$verbose_digits,
    verbose_show_ests = control.params$verbose_show_ests,
    tot_iter=iter
  )

  ## start sampling ####
  for (s in 2:iter) {

    ## continuous version of outcome (latent outcome under binomial probit model)
    if (est.h) {
      chain$ystar[s,] <- ystar.update(y = y, X = X, beta = chain$beta[s - 1,], h = chain$h[s - 1, ])
    } else {
      #Try without h, if it fails, use h
      res <- ystar.update.noh(y = y, X = X, beta = chain$beta[s - 1,], Vinv = Vcomps$Vinv, ystar = chain$ystar[s - 1, ])
      if(any(is.na(res))){
        if(all(chain$h[s-1,] == 0)){ #Make sure h is updated for the necessary iteration
          chain$h[s-1,] <- h.update(lambda = chain$lambda[s-1,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s-1], y = ycont, X = X, beta = chain$beta[s-1,], r = chain$r[s-1,], Z = Z, data.comps = data.comps)$hsamp
        }
        res <- ystar.update(y = y, X = X, beta = chain$beta[s - 1,], h = chain$h[s - 1, ])
      }
      chain$ystar[s,] <- res
    }
    ycont <- chain$ystar[s, ]

    ## generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)

    ## beta
    if (!missingX) {
      chain$beta[s,] <- beta.update(X = X, Vinv = Vcomps$Vinv, y = ycont, sigsq.eps = chain$sigsq.eps[s - 1])
    }

    ## lambda
    lambdaSim <- chain$lambda[s - 1,]
    for (comp in 1:data.comps$nlambda) {
      varcomps <- lambda.update(r = chain$r[s - 1,], delta = chain$delta[s - 1,], lambda = lambdaSim, whichcomp = comp, y = ycont, X = X, Z = Z, beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, data.comps = data.comps, control.params = control.params)
      lambdaSim <- varcomps$lambda
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.lambda[s,comp] <- varcomps$acc
      }
    }
    chain$lambda[s,] <- lambdaSim


    ## r
    rSim <- chain$r[s - 1,]
    comp <- which(!1:ncol(Z) %in% ztest)
    if (length(comp) != 0) {
      if (rmethod == "equal") { ## common r for those variables not being selected
        varcomps <- r.update(r = rSim, whichcomp = comp, delta = chain$delta[s - 1,],
                             lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,],
                             sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z,
                             data.comps = data.comps, control.params = control.params,
                             rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
                             rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
                             rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
                             rprop.logdens = rprop.logdens)
        rSim <- varcomps$r
        if (varcomps$acc) {
          Vcomps <- varcomps$Vcomps
          chain$acc.r[s, comp] <- varcomps$acc
        }
      } else if (rmethod == "varying") { ## allow a different r_m
        for (whichr in comp) {
          varcomps <- r.update(r = rSim, whichcomp = whichr, delta = chain$delta[s - 1,],
                               lambda = chain$lambda[s,], y = ycont, X = X, beta = chain$beta[s,],
                               sigsq.eps = chain$sigsq.eps[s], Vcomps = Vcomps, Z = Z,
                               data.comps = data.comps, control.params = control.params,
                               rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
                               rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
                               rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
                               rprop.logdens = rprop.logdens)
          rSim <- varcomps$r
          if (varcomps$acc) {
            Vcomps <- varcomps$Vcomps
            chain$acc.r[s, whichr] <- varcomps$acc
          }
        }
      }
    }
    ## for those variables being selected: joint posterior of (r,delta)
    if (varsel) {
      varcomps <- rdelta.update(r = rSim, delta = chain$delta[s - 1,],
                                lambda = chain$lambda[s,], y = ycont, X = X,
                                beta = chain$beta[s,], sigsq.eps = chain$sigsq.eps[s],
                                Vcomps = Vcomps, Z = Z, ztest = ztest,
                                data.comps = data.comps, control.params = control.params,
                                rprior.logdens = rprior.logdens, rprop.gen1 = rprop.gen1,
                                rprop.logdens1 = rprop.logdens1, rprop.gen2 = rprop.gen2,
                                rprop.logdens2 = rprop.logdens2, rprop.gen = rprop.gen,
                                rprop.logdens = rprop.logdens)
      chain$delta[s,] <- varcomps$delta
      rSim <- varcomps$r
      chain$move.type[s] <- varcomps$move.type
      if (varcomps$acc) {
        Vcomps <- varcomps$Vcomps
        chain$acc.rdelta[s] <- varcomps$acc
      }
    }
    chain$r[s,] <- rSim

    ###################################################
    ## generate posterior sample of h(z) from its posterior P(h | beta, sigsq.eps, lambda, r, y)

    if (est.h) {
      hcomps <- h.update(lambda = chain$lambda[s,], Vcomps = Vcomps, sigsq.eps = chain$sigsq.eps[s], y = ycont, X = X, beta = chain$beta[s,], r = chain$r[s,], Z = Z, data.comps = data.comps)
      chain$h.hat[s,] <- hcomps$hsamp
      if (!is.null(hcomps$hsamp.star)) { ## GPP
        Vcomps$hsamp.star <- hcomps$hsamp.star
      }
      rm(hcomps)
    }

    ###################################################
    ## generate posterior samples of h(Znew) from its posterior P(hnew | beta, sigsq.eps, lambda, r, y)

    if (!is.null(Znew)) {
      chain$hnew[s,] <- newh.update(Z = Z, Znew = Znew, Vcomps = Vcomps, lambda = chain$lambda[s,], sigsq.eps = chain$sigsq.eps[s], r = chain$r[s,], y = ycont, X = X, beta = chain$beta[s,], data.comps = data.comps)
    }

    ###################################################
    ## print details of the model fit so far
    print_diagnostics(verbose = verbose, opts = opts, curr_iter = s, tot_iter = iter, chain = chain, varsel = varsel, hier_varsel = hier_varsel, ztest = ztest, Z = Z, groups = groups)
  }

  #Remove r.params from control.params for later functions
  control.params$r.params <- NULL

  #Put data into chain return (so subsets are easily identified)
  chain$y <- y
  chain$X <- X
  chain$Z <- Z

  return(list(sampOutput = chain,
              starting.values = starting.values,
              control.params = control.params))
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
                            starting.values = list(rho = 1),
                            control.params,
                            nchains = 1,
                            iter = 1000,
                            warmup = 50,
                            ...){
  #####################
  #Set Starting Values
  #####################
  starting.values.default <- list(rho = 1)
  starting.values <- modifyList(starting.values.default, as.list(starting.values))


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

  #################################
  #Add Unused parameters to samples
  #################################
  samples$sigsq.eps <- rep(1, nchains*iter)
  samples$delta <- matrix(1, nchains*iter, ncol(Z))
  samples$r <- matrix(1, nchains*iter, ncol(Z))

  return(list(sampOutput = samples,
              starting.values = starting.values,
              control.params = control.params))
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
                              control.params = list(),
                              nchains = 1,
                              iter = 1000,
                              warmup = 50,
                              ...){
  #####################
  #Set Starting Values
  #####################
  starting.values.default <- list(rho = 2)
  starting.values <- modifyList(starting.values.default, as.list(starting.values))

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
  return(list(sampOutput = samples,
              starting.values = starting.values,
              control.params = control.params))
}
