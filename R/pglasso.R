#' Performs positive graphical lasso or standard graphical lasso by optimizing their dual problems.
#'
#' This function implements a simple block-coordinate descent algorithm to find the maximum of the regularized
#' Gaussiann log-likelihood. The  penalty  term involves the negative partial correlations.
#' @param S the sample covariance matrix
#' @param rho positive penalty (can be Inf)
#' @param L matrix of lower penalties (can be -Inf)
#' @param U matrix of upper penalties (can be Inf)
#' @param tol the convergence tolerance (default tol=1e-8)
#' @param pos.constr if TRUE (default) penalizes positive K_ij, if FALSE performs the standard dual graphical lasso.
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' print(TRUE)
#' 
##### Algorithm 3
pglasso <- function(S,rho=NULL, L=NULL,U=NULL,tol=1e-7,pos.constr=TRUE){
  d <- nrow(S)
  if (is.null(rho)==FALSE){
    if (pos.constr==FALSE){
      cat("** The algorithm maximizes the penalized log-likelihood function with the standard glasso penalty.\n")
      L <- -rho*(matrix(1,d,d)-diag(d))
      U <- rho*(matrix(1,d,d)-diag(d))
    } else {
      cat("** The algorithm maximizes the penalized log-likelihood function with the positive glasso penalty.\n")
      L <- matrix(0,d,d)
      U <- rho*(matrix(1,d,d)-diag(d))
    }
  } else {
    if (is.null(L)||is.null(U)){
      cat("Error: You need to specify  the  penalty parameter(s).\n")
      return()
    }
    cat("** The algorithm maximizes the penalized log-likelihood function with the general LU-penalty.\n")
  }
  # for computation of the dual gap we need to replace Inf with a large number etc
  # in all other parts of the procedure we use L and U normally
  L0 <- L
  L0[which(L0==-Inf)]  <- -1e+5
  U0 <- U
  U0[which(U0==Inf)]  <- 1e+5
  
  #compute the starting point
  if (min(eigen(S)$values)>0){
    cat("The algorithm starts at the sample covariance matrix.\n")
    Sig <- S
  } else{
    cat("S is not positive definite. Computing the starting point..")
    if (pos.constr==FALSE || (min(U+diag(d))>0 && max(L-diag(d))<0)){
      Z <- diag(diag(S))
    } else {
      if (min(U+diag(d))==0){
        cat("\n \n **Warning: This combination of L,U is not supported unless S is PD..\n")
        return()
      }
      Z <- Zmatrix(S)
    }
    cat("..")
    t <- 1
    while(!(min(L<=round(t*(Z-S),8)) && min(round(t*(Z-S),8)<=U))){
      t <- t/2
    }
    Sig <- (1-t)*S+t*Z
    Sig <- (1-t)*S+t*Z # this is the starting point
    cat(" DONE\n")
  }
  K <- solve(Sig)
  it <- 0
  cat("\n The algorithm will stop when the dual gap is below: ",tol,"\b.\n\n")
  cat("Iteration | Dual Gap\n")
  dualgap <- Inf
  while(dualgap > tol){
    for (j in 1:d){
      A <- rbind(diag(d-1),-diag(d-1))
      b <- c(S[j,-j]+L[j,-j],-(S[j,-j]+U[j,-j]))
      active <- which((b>-Inf & b<Inf)) 
      A <- A[active,]
      b <- b[active]
      y <- quadprog::solve.QP(Dmat=solve(Sig[-j,-j]),dvec=rep(0,d-1),Amat=t(A),bvec=b)$solution
      Sig[j,-j] <- Sig[-j,j] <- y
    }
    it <- it+1
    K <- solve(Sig)
    dualgap <- sum(diag(S%*%K))-d+sum(pmax(L0*K,U0*K)) 
    cat(it,"\t  | ",dualgap,"\n")
  }
  return(list(K=(K+t(K))/2,it=it))
}

