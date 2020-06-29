#' A coordinate descent algorithm to find the dual maximum likelihood estimator for associated Gaussian distributions
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param W inverse of the sample covariance matrix ???
#' @param edges the list of edges that should be updated
#' @param tol the convergence tolerance (default tol=1e-8)
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
LADualMLE <- function(W,edges=combn(1:nrow(W),2,simplify=FALSE), tol=1e-8){
  p <- ncol(W)
  row.names(W) <- colnames(W) <- as.character(1:p)
  Sig <- W
  Sig0 <- 2*Sig
  it <- 0
  if (p==2){
    return(list(Sig=solve(matrix(c(W[1,1],min(c(W[1,2],0)),min(c(W[1,2],0)),W[2,2]),2,2))))
  }
  while (sum(abs(Sig-Sig0))>tol){
    Sig0 <- Sig
    it <- it+1
    for (A in edges){
      i <- A[1]
      j <- A[2]
      B <- setdiff(1:p,A)
      L <- Sig[A,B]%*%solve(Sig[B,B])%*%Sig[B,A]
      if (L[1,2]>=W[i,j]/(W[i,i]*W[j,j]-W[i,j]^2)) Sig[A,A] <- solve(W[A,A])+L
      else {
        au <- (1+sqrt(1+4*L[1,2]^2*W[i,i]*W[j,j]))/2
        Sig[A,A] <- diag(c(L[1,1]+au/W[i,i],L[2,2]+au/W[j,j]))
      }
    }
    # for (i in 1:(p-1)){
    #   for (j in (i+1):p){
    #     A <- union(i,j)
    #     B <- setdiff(1:p,A)
    #     L <- Sig[A,B]%*%solve(Sig[B,B])%*%Sig[B,A]
    #     if (L[1,2]>=W[i,j]/(W[i,i]*W[j,j]-W[i,j]^2)) Sig[A,A] <- solve(W[A,A])+L
    #     else {
    #       au <- (1+sqrt(1+4*L[1,2]^2*W[i,i]*W[j,j]))/2
    #       Sig[A,A] <- diag(c(L[1,1]+au/W[i,i],L[2,2]+au/W[j,j]))
    #     }
    #   }
    # }
  }
  # output
  return(list(Sig=Sig,it))
}
