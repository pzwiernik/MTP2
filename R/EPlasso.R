#' Performs l1-regularized estimation of an inverse covariance matrix assuming nonnegativity of correlations
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param S the sample covariance matrix
#' @param lambda positive penalty
#' @param tol the convergence tolerance (default tol=1e-8)
#' @param pos.constr if TRUE (default) assumes nnonnegative correlations, if FALSE performs the standard dual graphical lasso.
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
##### Algorithm 3
EPlasso <- function(S,lambda,tol=1e-7,pos.constr=TRUE){
  d <- nrow(S)
  BM <- 1.0*(S>=0)
  CM <- 1.0*(S>=lambda)
  Sig <- S+lambda*diag(d) # starting point (the diagonal matches the optimum)
  K <- solve(Sig)
  tol <- 1e-6
  while(sum(diag(S%*%K))-d+lambda*sum(abs(K)) > tol){
    for (j in 1:d){
      if (pos.constr==TRUE){
        AA <- which(BM[j,]==0) # these are the entries with negative S[i,j]
        BB <- setdiff(setdiff(which(BM[j,]==1),j),(which(CM[j,]==1))) # small positive entries of S
        CC <- setdiff(which(CM[j,]==1),j) # large positive entries of S
        A <- rbind(diag(d)[AA,-j],diag(d)[BB,-j],-diag(d)[BB,-j],diag(d)[CC,-j],-diag(d)[CC,-j])
        b <- c(rep(0,length(AA)),S[j,BB],-lambda-S[j,BB],-lambda+S[j,CC],-lambda-S[j,CC])
        meq=length(AA)
      } else {
        A <- rbind(diag(d-1),-diag(d-1))
        b <- c(-lambda+S[j,-j],-lambda-S[j,-j])
        meq=0
      }
      y <- solve.QP(Dmat=solve(Sig[-j,-j]),dvec=rep(0,d-1),Amat=t(A),bvec=b,meq=meq)$solution
      Sig[j,-j] <- Sig[-j,j] <- y
    }
    K <- solve(Sig)
    print(sum(diag(S%*%K))-d+lambda*sum(abs(K)))
  }
  return(list(K=K,it=it))
}

