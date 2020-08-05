#' Performs Graphical Oriented LAZy Optimization by optimizing the dual problem.
#'
#' This function implements a simple block-coordinate descent algorithm to find the maximum of the regularized
#' Gaussiann log-likelihood  with  a an assymetric penalty of lasso type.
#' @param S the sample covariance matrix
#' @param rho positive penalty (can be Inf)
#' @param L matrix of lower penalties (can be -Inf)
#' @param U matrix of upper penalties (can be Inf)
#' @param tol the convergence tolerance (default tol=1e-8)
#' @param pos.constr if TRUE (default) penalizes positive K_ij, if FALSE performs the standard dual graphical lasso.
#' @param output if TRUE (default) the output will be printed.
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' print(TRUE)
#' 
##### Algorithm 3
golazo <- function(S,rho=NULL, L=NULL,U=NULL,tol=1e-7,pos.constr=TRUE,output=TRUE){
  d <- nrow(S)
  if (is.null(rho)==FALSE){
    if (pos.constr==FALSE){
      if (output==TRUE){
        cat("** The algorithm maximizes the penalized log-likelihood function with the standard glasso penalty.\n")
      }
      L <- -rho*(matrix(1,d,d)-diag(d))
      U <- rho*(matrix(1,d,d)-diag(d))
    } else {
      if (output==TRUE){
        cat("** The algorithm maximizes the penalized log-likelihood function with the positive glasso penalty.\n")
      }
      L <- matrix(0,d,d)
      U <- rho*(matrix(1,d,d)-diag(d))
    }
  } else {
    if (is.null(L)||is.null(U)){
      if (output==TRUE){
        cat("Error: You need to specify  the  penalty parameter(s).\n")
      }
      return()
    }
    if (output==TRUE){
      cat("** The algorithm maximizes the penalized log-likelihood function with the general LU-penalty.\n")
    }
  }
  # for computation of the dual gap we need to replace Inf with a large number etc
  # in all other parts of the procedure we use L and U normally
  L0 <- pmax(L,-sqrt(outer(diag(S+U),diag(S+U)))-S)
  diag(L0) <- diag(L)
  U0 <- pmin(U,sqrt(outer(diag(S+U),diag(S+U)))-S)
  diag(U0) <- diag(U)

  #compute the starting point
  if (min(eigen(S)$values)>1e-4){
    if (output==TRUE){
      cat("The algorithm starts at the sample covariance matrix.\n")
    }
    Sig <- S
  } else{
    if (output==TRUE){
      cat("S is not positive definite. Computing the starting point..")
    }
    if (pos.constr==FALSE || (min(U+diag(d))>0 && max(L-diag(d))<0)){
      Z <- diag(diag(S))
    } else {
      if (min(U+diag(d))==0){
        cat("\n \n **Warning: This combination of L,U is not supported unless S is PD..\n")
        return()
      }
      Z <- Zmatrix(S)
    }
    if (output==TRUE){
      cat("..")
    }
    t <- 1
    while(!(min(L<=round(t*(Z-S),8)) && min(round(t*(Z-S),8)<=U))){
      t <- t/2
    }
    Sig <- (1-t)*S+t*Z
    Sig <- (1-t)*S+t*Z # this is the starting point
    if (output==TRUE){
      cat(" DONE\n")
    }
  }
  K <- solve(Sig)
  it <- 0
  if (output==TRUE){
    cat("\n The algorithm will stop when the dual gap is below: ",tol,"\b.\n\n")
    cat("Iteration | Dual Gap\n")
  }
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
    #roundK <- K * (abs(K)>tol)
    dualgap <- sum(S*K)-d+sum(pmax(L0*(abs(K)>tol)*K,U0*(abs(K)>tol)*K)) 
    if (output==TRUE){
      cat(it,"\t  | ",dualgap,"\n")
    }
  }
  return(list(K=(K+t(K))/2,Sig=(Sig+t(Sig))/2,it=it))
}

