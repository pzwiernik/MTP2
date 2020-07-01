#' A coordinate descent algorithm to find the dual maximum likelihood estimator for locally associated Gaussian distributions
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param W inverse of the sample covariance matrix ???
#' @param tol the convergence tolerance (default tol=1e-8)
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
LADualMLE <- function(W, tol=1e-8){
  AM <- (abs(W)>1e-8)
  K <- W
  d <- nrow(K)
  it <- 0
  cat("** The algorithm maximizes the dual log-likelihood with edge positivity constraints.\n\n")
  cat("* To establish convergence we track the  KKT conditions. Dual feasibility holds at each step.\n")
  cat("  All edge covariances must be >=",-tol,"\b.    Complementary slackness in each entry <=",tol,"\b.\n")
  cat("Min. edge cov | Complementary slackness\n")
  Shat <- solve(K)
  zeros <- list()
  nonzeros <- list()
  for (j in 1:d){
    zeros[[j]] <- setdiff(which(abs(W[j,])<1e-8),j)
    nonzeros[[j]] <- setdiff(which(abs(W[j,])>=1e-8),j)
  }
  while (min(Shat*AM)< -tol || max(abs(Shat*(W-K)))>tol){
    for (j in 1:d){
      A <- rbind(diag(d)[zeros[[j]],-j],-diag(d)[nonzeros[[j]],-j])
      b <- c(W[j,zeros[[j]]],-W[j,nonzeros[[j]]])
      y <- solve.QP(Dmat=solve(K[-j,-j]),dvec=rep(0,d-1),Amat=t(A),bvec=b,meq=length(zeros[[j]]))$solution
      K[j,-j] <- K[-j,j] <- y
    }
    Shat  <- solve(K)
    cat(min(Shat*AM),"       |",max(abs(Shat*(W-K))),"\n")
  }
  return(list(Sig=(Shat+t(Shat))/2,it))
}
