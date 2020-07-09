#' The W matrix
#'
#' Computes the W matrix as introduced in Section 4.2 of Lauritzen, Uhler, Zwiernik (2017).
#' @param S a sample covariance matrix
#' @keywords shortest path, W matrix
#' @export
#' @examples
#' print(TRUE)
#' 
Wmatrix <- function(S){
  p <- nrow(S)
  R <- stats::cov2cor(S)
  # Compute the distances. Non-positive correlations correspond to very big distances.
  D <- matrix(0,p,p)
    for (i in 1:(p-1)){
      for (j in (i+1):p){
        if (R[i,j]>0) D[i,j] <- D[j,i] <- -log(R[i,j])
      }
    }
  ig <- igraph::graph.adjacency(D, mode="undirected", weighted=TRUE)
  #plot(ig, edge.label=round(E(ig)$weight, 3))
  W <- exp(-igraph::shortest.paths(ig))
  return(diag(sqrt(diag(S)))%*%W%*%diag(sqrt(diag(S))))
}


