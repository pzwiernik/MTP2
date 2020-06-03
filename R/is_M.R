#' Checks if a positive definite matrix is an M-matrix
#'
#' Checks if all entries of K are smaller than a small positive number tol 
#' @param K an inverse covariance matrix
#' @param tol numbers greater than tol are treated as positive
#' @keywords M-matrix
#' @export
#' @examples
#' 
#' 
is.M <- function(K,tol=1e-10){
  p <- nrow(K)
  return(prod(1*(K<= 1e-10)+diag(p)))
}