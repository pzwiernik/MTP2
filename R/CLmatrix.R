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
#' S <- matrix(c(11.376,0.745,8.446,2.066,7.686, -0.759,0.745, 33.069,0.672, 
#' 36.047,2.016, 32.059,8.446,0.672,8.937,0.315,6.864, -0.598,2.066, 36.047,
#' 0.315, 51.939,2.187, 41.587,7.686,2.016,6.864,2.187,7.645,0.384, -0.759, 
#' 32.059, -0.598, 41.587,0.384, 41.560),6,6)
#' C <- CLmatrix(S)
CLmatrix <- function(S,positive=TRUE){
  p <- nrow(S)
  R <- stats::cov2cor(S)
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

