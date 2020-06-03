#' A coordinate descent algorithm to find the maximum likelihood estimator for Markov positive distribution by local updates .
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param cliques the list ov vectors with indices of cliques
#' @param S the sample covariance matrix
#' @param n the sample size (default 1), relevant for testing 
#' @param tol the convergence tolerance (default tol=1e-8)
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
##### Algorithm 3
MPFit <- function(cliques,S,n=1,tol=1e-8){
  d <- nrow(S)
  cliques <- lapply(cliques,as.character)
  row.names(S) <- colnames(S) <- as.character(1:d)
  K <- gRim::ggmfit(S,n.obs=n,glist= cliques)$K
  V <- 1:d
  K0 <- 2*K
  it <- 0
  while (sum(abs(K-K0))>tol){
    it <- it+1
    K0 <- K
    for (C in cliques){
      LC <- as.matrix(coord2K(S[C,C])$K)
      K[C,C] <- LC+K[C,setdiff(V,C)]%*%solve(K[setdiff(V,C),setdiff(V,C)])%*%K[setdiff(V,C),C]
    }
  }
  return(list(K=K,it=it))
}
