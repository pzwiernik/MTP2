#' Checks if a positive definite matrix is an M-matrix
#'
#' Checks if all entries of K are smaller than a small positive number tol 
#' @param S the sample covariance matrix
#' @param K an M-matrix
#' @param tol numbers greater than tol are treated as positive
#' @keywords 
#' @export
#' @examples
#' 
#' 
is.M <- function(K,tol=1e-10){
  p <- nrow(K)
  return(prod(1*(K<= 1e-10)+diag(p)))
}