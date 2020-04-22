#' The (positive) Chow-Liu matrix
#'
#' Computes the Chow-Liu estimate or the positive Chow-Liu estimate in which 
#' case only nonnnegative entries of S are taken into account. The Chow-Liu estimate 
#' corresponds to the maximum likelihood over all Gaussian tree models.
#' @param S a sample covariance matrix
#' @param positive TRUE
#' @keywords Chow-Liu algorithm, minimum cost spanning tree
#' @export
#' @examples
#' data(carcass)
#' S <- cov(carcass)[1:6,1:6]
#' C <- CLmatrix(S)
#' graphK(solve(C))
CLmatrix <- function(S,positive=TRUE){
  p <- nrow(S)
  R <- cov2cor(S)
  T <- MWSF(S,positive)  #plot(ig, edge.label=round(E(ig)$weight, 3))
  if (positive==TRUE) {W <- exp(-igraph::shortest.paths(T))}
  else {
    W <- diag(p)
    for (i in 1:p){
      for (j in setdiff(1:p,i)) {
        indcs <- as.vector(igraph::all_simple_paths(T,i,j)[[1]])
        l <- length(indcs)-1
        W[i,j] <- 1
        for (k in 1:l) {W[i,j] <- W[i,j]*R[indcs[k],indcs[k+1]]}
      }
    }
  }
return(diag(sqrt(diag(S)))%*%W%*%diag(sqrt(diag(S))))
  }

