#' The Chow-Liu sign adjustment
#'
#' Computes the Chow-Liu estimate or the positive Chow-Liu estimate in which 
#' case only nonnnegative entries of S are taken into account. The Chow-Liu estimate 
#' corresponds to the maximum likelihood over all Gaussian tree models.
#' @param S a sample covariance matrix
#' @keywords Chow-Liu algorithm, minimum cost spanning tree
#' @export
#' @examples 
#' print(TRUE)
#' 
#' 
#' 
CLsigns <- function(S){
  p <- nrow(S)
  T <- MWSF(S,FALSE)
  SS <- rep(0,p) # our D matrix
  SS[1] <- 1 # we will root the tree in 1
  for (i in 2:p){
      indcs <- as.vector(igraph::all_simple_paths(T,1,i)[[1]])
      l <- length(indcs)
      for (k in 2:l) {SS[indcs[k]] <- sign(SS[indcs[k-1]]*S[indcs[k-1],indcs[k]])}
  }
  return(SS)
  }
  