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
#' library(gRbase)
#' data(carcass)
#' S <- cov(carcass)[1:6,1:6]
#' C <- CLmatrix(S)
#' graphK(solve(C))
MWSF <- function(S,positive=TRUE){
  p <- nrow(S)
  R <- cov2cor(S)
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

