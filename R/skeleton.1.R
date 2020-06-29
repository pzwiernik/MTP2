#' Compute the 1-skeleton of a simplicial complex
#'
#' This function implements a simple coordinate descent algorithm to find the maximum likelihood estimator over
#' Gaussian MTP2 distributions. For details see Lauritzen, Uhler, Zwiernik (2017).
#' @param simplices the list ov maximal simplices
#' @return the optimal value of the concentration matrix
#' @return the number of iterations the algorithm needed to converge
#' @return the corresponding value of the log-likelihood
#' @keywords coordinate descent, concentration matrix.
#' @export
#' @examples
#' 
#' 
skeleton.1 <- function(simplices){
  skel <- list()
  for (C in simplices){
    skel <- append(skel,combn(C,2,simplify=FALSE))
  }
  return(unique(skel))
}
