#' Coordinate descent algorithm for the MLE.
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. The algorithm performs local updates of the covariance matrix $Sigma$.
#' For details see Lauritzen, Uhler, Zwiernik (2017); look for Algorithm 4.
#' @param S the sample covariance matrix
#' @param n the sample size (default 1)
#' @param Sig0 starting point of the algorithm
#' @param tol convergence criterion sum(abs(Sig-Sig0))<tol
#' @keywords coordinate descent, covariance matrix.
#' @export
#' @examples
#' library(gRbase)
#' data("carcass")
#' S <- cov(carcass)[1:6,1:6]
#' coord2Sig(S)
#'


##### Algorithm 4
coord2Sig <- function(S,n=1,Sig0=S,tol=1e-10,graph=FALSE){
  p <- nrow(S)
  Sig <- Sig0
  # start with G, the complete graph
  Sig0 <- 2*Sig
  it <- 0
#  R <- cov2cor(S)
#  Z <- Zmatrix(S)
  while (sum(abs(Sig-Sig0))>tol){
    Sig0 <- Sig
    it <- it+1
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        A <- union(i,j)
        B <- setdiff(1:p,A)
        L <- Sig[A,B]%*%solve(Sig[B,B])%*%Sig[B,A]
        #if (max(c(L[1,2],S[i,j]))-Sig[i,j]>0.001) print(max(c(L[1,2],S[i,j]))-Sig[i,j])
        Sig[i,j] <- Sig[j,i] <- max(c(L[1,2],S[i,j]))
       # if (is.M(solve(Sig))==1) print(it)
#        print(R[i,j]-1e-6<=Sig[i,j] && Sig[i,j]<=Z[i,j]+1e-6)
      }
    }
  }
  if (graph==TRUE) {graphK(solve(cov2cor(Sig)))}
  # output
  LR <- n*(-log(det(S))+log(det(Sig)))
  dimnames(Sig) <- dimnames(S)
  if (n==1) {
    return(list(Sig,it,(n/2)*(-log(det(Sig))-p),print("Sample size not provided")))
  }
  else return(list(Sig,it,(n/2)*(-log(det(Sig))-p),mtpLR(solve(Sig),LR)))
}
