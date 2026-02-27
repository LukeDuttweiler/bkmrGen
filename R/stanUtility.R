prepForStan <- function(y,
                        Z,
                        X,
                        starting.values,
                        M){
  dt <- starting.values
  dt$y <- y
  dt$Z <- Z
  dt$X <- X
  dt$p <- ncol(Z)
  dt$d <- ncol(X)
  dt$N <- length(y)
  dt$M <- M

  return(dt)
}
