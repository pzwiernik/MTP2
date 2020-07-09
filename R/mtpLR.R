#' The likelihood ratio test for MTP2
#'
#' Internal funciton that computes the likelihood ratio statistics to test the MTP constraints. It also returns the degrees 
#' of freedom and the corresponding p-value.
#' @param K an M-matrix that maximizes the likelihood
#' @param LR the value of the likelihood ratio statistics with respect to the unconstrained model
#' @param tol numbers greater than -tol in K are treated as zero (default 1e-8)
#' @keywords xxx
#' @export
#' @examples
#' print(TRUE)
mtpLR <- function(K,LR,tol=1e-10){
  p <- nrow(K)
  G <- igraph::graph_from_adjacency_matrix(1*(K< -tol),mode="undirected")
  df <- p*(p-1)/2-nrow(igraph::get.edges(G,igraph::E(G)))
  return(list(LR,df,stats::pchisq(LR,df)))
}