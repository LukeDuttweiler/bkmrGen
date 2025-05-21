#' Generate Trace Plots for Parameters from a BKMR Fit
#'
#' Produces trace plots for posterior samples of parameters from a fitted Bayesian Kernel Machine Regression (BKMR) model.
#' This function is designed to work with objects created using the \code{kmbayes()} function.
#'
#' @param bkmrFit An object of class \code{bkmrfit}, typically the result of the \code{kmbayes()} function containing MCMC samples.
#' @param par Character string specifying the parameter to plot the trace for.
#'   Valid options are \code{"beta"}, \code{"lambda"}, \code{"sigsq.eps"}, \code{"h.hat"}, \code{"r"}, or \code{"delta"}.
#' @param comp Integer(s) or \code{"all"} indicating which component(s) of the parameter to plot. If \code{"all"}, trace plots are generated for all components.
#' @param subset An integer indicating which posterior sample subset to use. Defaults to 1.
#'
#' @return A list of ggplot2 objects (or a single plot) containing the trace plot(s) for the specified parameter/component.
#'
#' @details This function extracts posterior samples for the chosen parameter and component(s), splits them by MCMC chain,
#' and generates trace plots using the helper function \code{genDiagnostic()}.
#' When multiple components are requested, the samples are reshaped into the format expected by \code{genDiagnostic()}.
#'
#' @examples
#' \dontrun{
#' fit <- kmbayes(...)  # assuming the fit has been run
#' TracePlot(fit, par = "beta", comp = 'all')         # Trace plot for all beta coefficiens
#' TracePlot(fit, par = "lambda", comp = 1)   # Trace plots for all lambda components
#' }
#'
TracePlot <- function(bkmrFit,
                      par = 'beta',
                      comp = 1,
                      subset = 1){


  #Ensure par is an actual option
  if(!(par %in% c('beta', 'lambda', 'sigsq.eps', 'h.hat', 'r', 'delta', 'rho'))){
    stop("par must be one of 'beta', 'lambda', 'sigsq.eps', 'h.hat', 'r', or 'delta'")
  }

  if(!is.numeric(comp)){
    if(comp == 'all'){
      comp <- 1:ncol(bkmrFit$subSamples[[subset]][[par]])
    }
  }

  #Retrieve desired parameters and components, separate into chains
  nchains <- bkmrFit$nchains
  iter <- bkmrFit$iter
  chainInd <- rep(1:nchains, each = iter)
  samps <- bkmrFit$subSamples[[subset]][[par]][,comp, drop = F]
  samps <- lapply(1:nchains, function(ch){
    return(samps[ch == chainInd,])
  })

  #Produce diagnostics with genMCMCDiag, different functionality depending on components
  if(length(comp) == 1){
    tp <- genMCMCDiag::genDiagnostic(samps, diagnostics = 'traceplot')
  }else{
    #Break samps into list for genDiagnostic structure
    samps <- lapply(samps, function(ch){
      return(as.list(as.data.frame(t(ch))))
    })

    tp <- genMCMCDiag::genDiagnostic(samps, diagnostic = 'traceplot',
                                     proximityMap = 'lanfear', distance = genMCMCDiag::eucDist,
                                     reference = 0)
  }

  return(tp)
}
