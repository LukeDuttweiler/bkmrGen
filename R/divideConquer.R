#' Get Univariate Wasserstein Posterior
#'
#' Takes list of samples from univariate sub-sample posteriors and returns samples from the Wasserstein Posterior. Utilizes Algorithm 1 fromSrivastava et al. (2018), assuming theta_ij is univariate and f is the identity function.
#'
#' @param sampList List of samples from univariate posterior. Each sample must be a numeric vector.
#' @param numAtoms Number of atoms desired for calculating WASP
#' @param solver Linear-program solver to use. Solvers available are provided by the R package \code{\link{ROI}}.
#' @param n_samps Number of samples to return from WASP. Defaults as length of first sub-sample provided in sampList.
#'
#' @return n_samps samples from WASP
#'
wasp_univariate <- function(sampList, numAtoms = 100, solver = 'lpsolve',
                            n_samps = length(sampList[[1]])){
  #IF only one sample is provided, return that sample
  if(length(sampList) == 1){
    return(sampList[[1]])
  }

  #Algorithm 1 from Srivastava et al. (2018), assuming \theta_ij is univariate and
  #f is the identity function

  #STEP 1
  phiMin <- min(sapply(sampList, min))
  phiMax <- max(sapply(sampList, max))

  #STEP 2
  #g <- ceiling((phiMax - phiMin)/meshSize)

  #STEP 3
  atomMat <- seq(phiMin, phiMax, length.out =numAtoms)

  #STEP 4
  Djs <- lapply(sampList, function(samp){
    return(outer(atomMat, samp, FUN = '-')^2)
  })

  #STEP 5
  medPost <- wasp_linProg(DList = Djs, atomMat = atomMat, solver = solver)
  medPost <- pmax(medPost, rep(0, length(medPost)))

  #SAMPLE FROM ATOMS BASED ON SOLUTION
  wasp_sample <- sample(atomMat, length(sampList[[1]]), replace = TRUE, prob = medPost)

  #GET KERNEL DENSITY ESTIMATE (Not sure we need this?)
  #wasp_dens <- KernSmooth::bkde(wasp_sample, bandwidth = max(diff(atomMat)[1], 1),
  #                              range = range(atomMat))

  #RETURN POSTERIOR RESULTS
  return(wasp_sample)
}



#' Run Linear Program for calculating WASP
#'
#' Linear program for calculating Wasserstein Posterior (WASP). See \code{\link{wasp_univariate}}.
#'
#' @param DList List of D matrices from WASP algorithm
#' @param atomMat Matrix of atoms from WASP algorithm
#' @param solver Linear program solver to use, uses ROI package.
#'
#' @return Solution to WASP linear program
#'
wasp_linProg <- function(DList, atomMat, solver = 'lpsolve'){
  #constants
  K <- length(DList)
  Ni <- ncol(DList[[1]])
  N <- length(atomMat)
  In <- diag(N)
  En <- matrix(1, nrow = 1, ncol = N)

  #Build matrix A0, encodes constraint T_j^T1_g = (1_sj)/sj
  #A0 is block diagonal, each block corresponds to one of the T_j matrices
  sjFrac <- t(rep(1, Ni)) %x% In/Ni
  A0 <- Matrix::bdiag(replicate(K, sjFrac, simplify = FALSE))

  #Build matrix A00, encodes constraint T_j1_sj = a
  A00 <- rep(1,K) %x% -In

  #Matrix A1 captures all of these constraints
  A1 <- Matrix::drop0(cbind(A00, A0))

  #Matrix A2 captures all simplex constraints
  A2 <- Matrix::bdiag(replicate(Ni*K + 1, En, simplify = FALSE))

  #Full constraint matrix
  A <- rbind(A1, A2)

  #Change to slam 'simple_triplet_matrix' to match ROI input requirements
  #NOTE: This avoids converting back to a full matrix, so saves a TON of memory
  A <- Matrix_to_slam(A)

  #Right hand side of constraints, generated as b
  b <- c(rep(0, K*N), rep(1, Ni*K + 1))

  #Generate cost vector (costVec), first set of 0s is for cost of a values
  #Second set is cost of elements of T_ij
  costCell <- lapply(DList, function(x){x/ncol(x)})
  costVec <- c(rep(0, N), do.call(c, costCell))

  #Now prep for LP Solver
  #Objective function
  obj <- ROI::L_objective(costVec)

  #Constraints
  const <- ROI::L_constraint(L = A,
                             rhs = b,
                             dir = rep('==', length(b)))

  #LP Model
  lp <- ROI::OP(obj, const)

  #Solution
  solution <- ROI::ROI_solve(lp, solver = solver)

  return(ROI::solution(solution)[1:length(atomMat)])
}
