#' A coordinate descent algorithm to find the maximum likelihood estimator by local updates of the concntration matrix K.
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param S the sample covariance matrix
#' @param n the sample size (default 1), relevant for testing 
#' @param K0 the starting point (default K0=solve(diag(diag(S))))
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
coord2K <- function(S,n=1,K0=solve(diag(diag(S))),tol=1e-8){
p <- ncol(S)
K <- K0
# start with G, the complete graph
K0 <- 2*K
it <- 0

while (sum(abs(K-K0))>tol){
K0 <- K
it <- it+1
for (i in 1:(p-1)){
  for (j in (i+1):p){
    A <- union(i,j)
    B <- setdiff(1:p,A)
    L <- K[A,B]%*%solve(K[B,B])%*%K[B,A]
    if (L[1,2]<=S[i,j]/(S[i,i]*S[j,j]-S[i,j]^2)) K[A,A] <- solve(S[A,A])+L
    else {
      au <- (1+sqrt(1+4*L[1,2]^2*S[i,i]*S[j,j]))/2
      K[A,A] <- diag(c(L[1,1]+au/S[i,i],L[2,2]+au/S[j,j]))
     }
  }
}
}
# output
LR <- n*(-log(det(S))-log(det(K)))
dimnames(K) <- dimnames(S)
if (n==1) {
  return(list(K,it,(n/2)*(log(det(K))-sum(diag(S%*%K))),print("Sample size not provided")))
}
else return(list(K,it,(n/2)*(log(det(K))-sum(diag(S%*%K))),mtpLR(K,LR)))
}
