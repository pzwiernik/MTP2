#' A coordinate descent algorithm to find the maximum likelihood estimator for Markov positive distribution by local updates .
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param cliques the list ov vectors with indices of cliques
#' @param S the sample covariance matrix
#' @param n the sample size (default 1), relevant for testing 
#' @param tol the convergence tolerance (default tol=1e-8)
#' @param skelstart if TRUE first a 1-skeleton model is fitted and used as the starting point (default FALSE) 
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
##### Algorithm 3
MPFit <- function(cliques,S,n=1,tol=1e-7,skelstart=FALSE){
  d <- nrow(S)
  cliques <- lapply(cliques,as.character)
  Z <- Zmatrix(S)
  row.names(S) <- colnames(S) <- row.names(Z) <- colnames(Z) <- as.character(1:d)
  K <- gRim::ggmfit(Z,n.obs=n,glist= cliques)$K
  V <- 1:d
  K0 <- 2*K
  it <- 0
  kl <- c()
  KD <- solve(Z)
  if (skelstart == TRUE){
    cliques0 <- skeleton.1(cliques)
    while (MP.KKT.check(cliques0,S,K,tol)==0){
      it <- it+1
      K0 <- K
      for (C in cliques0){
        LC <- as.matrix(coord2K(S[C,C])$K)
        B <- setdiff(V,C)
        K[C,C] <- (LC+K[C,B]%*%solve(K[B,B])%*%K[B,C])
        kl <- c(kl,KL(S,solve(KD))-KL(S,solve(K))-KL(solve(K),solve(KD)))
        KD <- K
      }
    }
    print(K)
  }
  while (MP.KKT.check(cliques,S,K,tol)==0){
    it <- it+1
    K0 <- K
    for (C in cliques){
      LC <- as.matrix(coord2K(S[C,C])$K)
      B <- setdiff(V,C)
      K[C,C] <- (LC+K[C,B]%*%solve(K[B,B])%*%K[B,C])
      kl <- c(kl,KL(S,solve(KD))-KL(S,solve(K))-KL(solve(K),solve(KD)))
      #print(KL(S,solve(K1))-KL(S,solve(K))-KL(solve(K),solve(K1)))
      print(K)
      print(solve(solve(K)[c(1,2,3),c(1,2,3)]))
      print(solve(solve(K)[c(2,3,4),c(2,3,4)]))
      KD <- K
    }
  }
  return(list(K=K,it=it,kl=kl))
}

MP.KKT.check <- function(cliques,S,K,tol=1e-7){
  Sig <- solve(K)
  check <- 1
  for (C in cliques){
    LC <- solve(Sig[C,C])
    check <- check*is.M(LC,tol)
    check <- check*min(Sig[C,C]-S[C,C]+tol>=0)
    check <- check*min(abs(LC*(Sig[C,C]-S[C,C]))<tol)
  }
  return(check)
}
