#' The graph of an M-matrix
#'
#' Plots the graph representing the support of an M-matrix. 
#' Small negative numbers are treated as zero.
#' @param S the sample covariance matrix
#' @param K an M-matrix
#' @param tol numbers greater than -tol are treated as zero
#' @keywords xxx
#' @export
#' @examples
#' 
#' 
is.M <- function(K,tol=1e-10){
  p <- nrow(K)
  return(prod(1*(K<= 1e-10)+diag(p)))
}