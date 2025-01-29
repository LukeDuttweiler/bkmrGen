wasp_univariate(samp, indx, meshSize = 100){
  #Algorithm 1 from Srivastava et al. (2018), assuming \theta_ij is univariate and
  #f is the identity function

  #STEP 1
  phiMin <- min(samp)
  phiMax <- max(samp)

  #STEP 2
  g <- ceiling((phiMax - phiMin)/meshSize)

  #STEP 3
  atomMat <- matrix(phiMin + ((1:g)/g)*(phiMax - phiMin), ncol = 1)

  #STEP 4
  Djs <- lapply(1:max(indx), function(i){
    si <- samp[indx == i]
    return(outer(atomMat, si, FUN = '-')^2)
  })
}
