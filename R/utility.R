#' Convert 'Matrix' object to 'slam' object
#'
#' Converts a 'Matrix' object to a 'slam' object without reverting to 'matrix' type. This saves a ton of memory and time.
#'
#' @param Mat Object from package 'Matrix'
#'
#' @return Matrix with type 'slam'
Matrix_to_slam <- function(Mat){
  p <- rep(seq_along(diff(Mat@p)), diff(Mat@p))
  Mat <- slam::simple_triplet_matrix(i = (Mat@i)+1, j = p, v = Mat@x,
                                   nrow = Mat@Dim[1], ncol = Mat@Dim[2])
  return(Mat)
}

#' Make K_{Z,r} matrix (Exponential kernel)
#'
#' See Bobb et al. (2015), supplement page 3 for details.
#'
#' @param r Vector of augmented variables in kernel matrix controlling smoothness of h()
#' @param Z1 Exposure design matrix
#' @param Z2 New exposure design matrix (optional)
#'
#' @return K_{Z,r} matrix
makeKpart <- function(r, Z1, Z2 = NULL) {
  Z1r <- t(t(Z1) * c(sqrt(r)))
  if (is.null(Z2)) {
    Kpart <- fields::rdist(Z1r)^2
  } else {
    Z2r <- t(t(Z2) * c(sqrt(r)))
    Kpart <- fields::rdist(Z1r, Z2r)^2
  }
  Kpart
}


#' Make necessary components from V_{lambda, Z, r}
#'
#' See Bobb et al. (2015), supplement page 4 for details.
#'
#' @param r Vector of augmented variables in kernel matrix controlling smoothness of h()
#' @param lambda Variance parameter == tau*sigma^{-2}
#' @param Z Exposure design matrix
#' @param data.comps List of control parameters to change behavior. Internally created by \code{\link{kmbayes}}
#'
#' @return List containing V_{lambda, Z, r}^{-1} and log(determinant(V)) and additional values if selecting knots for kernel
makeVcomps <- function(r, lambda, Z, data.comps) {
  if (is.null(data.comps$knots)) {
    Kpart <- makeKpart(r, Z)
    V <- diag(1, nrow(Z), nrow(Z)) + lambda[1]*exp(-Kpart)
    if (data.comps$nlambda == 2) {
      V <- V + lambda[2]*data.comps$crossTT
    }
    cholV <- chol(V)
    Vinv <- chol2inv(cholV)
    logdetVinv <- -2*sum(log(diag(cholV)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv)
  } else {## predictive process approach
    ## note: currently does not work with random intercept model
    nugget <- 0.001
    n0 <- nrow(Z)
    n1 <- nrow(data.comps$knots)
    nall <- n0 + n1
    K1 <- exp(-makeKpart(r, data.comps$knots))
    K10 <- exp(-makeKpart(r, data.comps$knots, Z))
    Q <- K1 + diag(nugget, n1, n1)
    R <- Q + lambda[1]*tcrossprod(K10)
    cholQ <- chol(Q)
    cholR <- chol(R)
    Qinv <- chol2inv(cholQ)
    Rinv <- chol2inv(cholR)
    Vinv <- diag(1, n0, n0) - lambda[1]*t(K10) %*% Rinv %*% K10
    logdetVinv <- 2*sum(log(diag(cholQ))) - 2*sum(log(diag(cholR)))
    Vcomps <- list(Vinv = Vinv, logdetVinv = logdetVinv, cholR = cholR, Q = Q, K10 = K10, Qinv = Qinv, Rinv = Rinv)
  }
  Vcomps
}

#' Create Values needed in makeVcomps, depending on if using random intercept model
#'
#' @param id Vector of id values (for random intercept model), default in kmbayes is NULL
#' @param knots Vector of knots for Gaussian Process modification of V
#'
#' @return List of options used by makeVcomps to estiamte V
createDataComps <- function(id, knots){
  if (!is.null(id)) { ## for random intercept model
    randint <- TRUE
    id <- as.numeric(as.factor(id))
    nid <- length(unique(id))
    nlambda <- 2

    ## matrix that multiplies the random intercept vector
    TT <- matrix(0, length(id), nid)
    for (i in 1:nid) {
      TT[which(id == i), i] <- 1
    }
    crossTT <- tcrossprod(TT)
    rm(TT, nid)
  } else {
    randint <- FALSE
    nlambda <- 1
    crossTT <- 0
  }
  data.comps <- list(randint = randint, nlambda = nlambda, crossTT = crossTT)
  if (!is.null(knots)) data.comps$knots <- knots

  return(data.comps)
}


#' Control bkmr printed output frequency and detail
#'
#' @param verbose_freq Frequency of output
#' @param verbose_show_ests Logical, show estimates or not
#' @param verbose_digits number of digits to print
#' @param tot_iter Number of iterations in function
#'
#' @return Returns list of options
set_verbose_opts <- function(verbose_freq = NULL, verbose_show_ests = NULL, verbose_digits = NULL,
                             tot_iter) {
  if (is.null(verbose_freq)) verbose_freq <- 10
  if (is.null(verbose_digits)) verbose_digits <- 5
  if (is.null(verbose_show_ests)) verbose_show_ests <- FALSE

  all_iter <- 100*(1:tot_iter)/tot_iter
  sel_iter <- seq(verbose_freq, 100, by = verbose_freq)
  print_iter <- sapply(sel_iter, function(x) min(which(all_iter >= x)))

  opts <- list(
    verbose_freq = verbose_freq,
    verbose_digits = verbose_digits,
    verbose_show_ests = verbose_show_ests,
    print_iter = print_iter
  )
  opts
}

loadROISolver <- function(solver){
  pckgName <- paste0('ROI.plugin.', solver)

  #Does the package exist?
  pckgExist <- pckgName %in% rownames(installed.packages())

  #If not, attempt to install
  if(!pckgExist){
    install.packages(pckgName)
  }

  #Then load
  eval(parse(text = paste0('library(', pckgName, ')')))

  return(invisible(NULL))
}
