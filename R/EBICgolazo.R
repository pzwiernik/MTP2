#' Finds  the optimal penalty parameter for Graphical Oriented Lasso using the EBIC criterion.
#'
#' This function computes the EBIC criterion for a seried of rho parameters. The penalty matrices are $rho L$, $rho U$, 
#' where L and U are fixed in  advance. 
#' @param S the sample covariance matrix
#' @param n the sample size
#' @param L matrix of lower penalties (can be -Inf)
#' @param U matrix of upper penalties (can be Inf)
#' @param tol the convergence tolerance (default tol=1e-8)
#' @param edge.tol which entries of K are treated as zero
#' @param gamma the EBIC parameter between 0 and 1
#' @param rhomin minimal rho to check
#' @param rhomax maximal rho to check
#' @param nrhos positive penalty (can be Inf)
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
EBICgolazo <- function(S,n=NULL,L= NULL,U=NULL,tol=1e-6,edge.tol=1e-4,gamma=0.5,rhomin=0.01,rhomax=1,nrhos=50,output=TRUE){
  # rhomax=1 makes a lot of sense if S is a correlation matrix
  d <- nrow(S)
  if (is.null(L)||is.null(U)){
      if (output==TRUE){
        cat("The value of L and U set to the default value.")
      }
     L <- matrix(0,d,d)
     U <- matrix(1,d,d)-diag(d)
  }
  rhos <- seq(from=rhomin,to=rhomax,length.out=nrhos)
  ebic <- rhos
  ebic.gamma <- gamma
  for (rho in rhos){
    LL <- rho*L
    UU <- rho*U
    res <- golazo(S,L=LL,U=UU,output=FALSE)
    # compute EBIC
    K <- res$K
    KR <- cov2cor(K) #to make edge count independend of scalings
    nedg <- length(which(abs(KR[upper.tri(abs(KR),diag=FALSE)])> edge.tol))
    ebic[which(rhos==rho)] <- -(n)*(log(det(K))-sum(S*K))+nedg * (log(n) + 4 * ebic.gamma * log(d))
  }
  return(list(rhos=rhos,ebic=ebic,opt.rho=rhos[min(which(ebic==min(ebic)))]))
}

