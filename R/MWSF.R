#' The (positive) maximum weight spanning forest
#'
#' Computes the Chow-Liu estimate or the positive Chow-Liu estimate in which 
#' case only nonnnegative entries of S are taken into account. The Chow-Liu estimate 
#' corresponds to the maximum likelihood over all Gaussian tree models.
#' @param S a sample covariance matrix
#' @param positive TRUE
#' @keywords Chow-Liu algorithm, minimum cost spanning tree
#' @export
#' @examples
#' S <- matrix(c(11.376,0.745,8.446,2.066,7.686, -0.759,0.745, 33.069,
#' 0.672, 36.047,2.016, 32.059,8.446,0.672,8.937,0.315,6.864, -0.598,
#' 2.066, 36.047,0.315, 51.939,2.187, 41.587,7.686,2.016,6.864,2.187,
#' 7.645,0.384, -0.759, 32.059, -0.598, 41.587,0.384, 41.560),6,6)
#' MWSF(S)
MWSF <- function(S,positive=TRUE){
  p <- nrow(S)
  R <- stats::cov2cor(S)
  # Compute the distances. Non-positive correlations correspond to very big distances.
  D <- matrix(0,p,p)
  if (positive==TRUE){
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        if (R[i,j]>0) D[i,j] <- D[j,i] <- -log(R[i,j])
      }
    }}
  else {  for (i in 1:(p-1)){
    for (j in (i+1):p){
      D[i,j] <- D[j,i] <- -log(abs(R[i,j]))
    }
  }}
  ig <- igraph::graph.adjacency(D, mode="undirected", weighted=TRUE)
  return(igraph::minimum.spanning.tree(ig))
}

