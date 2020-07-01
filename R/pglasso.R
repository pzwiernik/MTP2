#' Performs positive graphical lasso or standard graphical lasso by optimizing their dual problems.
#'
#' This function implements a simple block-coordinate descent algorithm to find the maximum of the regularized
#' Gaussiann log-likelihood. The  penalty  term involves the negative partial correlations.
#' @param S the sample covariance matrix
#' @param rho positive penalty (can be Inf)
#' @param tol the convergence tolerance (default tol=1e-8)
#' @param pos.constr if TRUE (default) penalizes positive K_ij, if FALSE performs the standard dual graphical lasso.
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
##### Algorithm 3
pglasso <- function(S,rho,tol=1e-7,pos.constr=TRUE){
  d <- nrow(S)
  if (pos.constr==FALSE){
    cat("** The algorithm maximizes the penalized log-likelihood function with the standard glasso penalty.\n")
  } else {
    cat("** The algorithm maximizes the penalized log-likelihood function with the positive glasso penalty.\n")
  }
  #compute the starting point
  cat("Computing the starting point..")
  if (pos.constr==FALSE){
    Z <- diag(diag(S))
} else {
    Z <- Zmatrix(S)
}
  cat("..")
  t <- 1
  mm <- max(abs(Z-S))
  if (mm>0 && rho<Inf){
    t <- min(1,rho/max(abs(Z-S)))
  }
  Sig <- (1-t)*S+t*Z # this is the starting point
  cat(" DONE\n")
  
  K <- solve(Sig)
  it <- 0
  cat("The algorithm stops when the dual gap is below: ",tol,"\b.\n\n")
  cat("Iteration | Dual Gap\n")
  dualgap <- Inf
  while(dualgap > tol || it==0){
    for (j in 1:d){
      A <- rbind(diag(d-1),-diag(d-1))
      if (pos.constr==TRUE){
        b <- c(S[j,-j],-rho-S[j,-j])
      } else{
        b <- c(S[j,-j]-rho,-rho-S[j,-j])
      }
      y <- solve.QP(Dmat=solve(Sig[-j,-j]),dvec=rep(0,d-1),Amat=t(A),bvec=b)$solution
      Sig[j,-j] <- Sig[-j,j] <- y
    }
    it <- it+1
    K <- solve(Sig)
    if (pos.constr==TRUE){
      dualgap <- sum(diag(S%*%K))-d+rho*sum((K-diag(diag(K)))*(K>0)) 
    } else{
      dualgap <- sum(diag(S%*%K))-d+rho*sum(abs(K-diag(diag(K)))) 
    }
    cat(it,"\t  | ",dualgap,"\n")
  }
  return(list(K=(K+t(K))/2,it=it))
}

