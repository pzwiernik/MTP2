#' Locally positive simplicial complex of a covariance matrix.
#' This function gives the set of all subblocks C for which Sigma[C,C] is an inverse M-matrix
#' @param S the sample covariance matrix
#' @param tol the convergence tolerance (default tol=1e-8)
#' @return the simplicial complex
#' @keywords inverse M-subblocks, simplicial complex, local positivity
#' @export
#' @examples
#' print(TRUE)
#' 
LP.simplicial.complex <- function(S,tol=1e-7){
  d <- nrow(S)
  V <- 1:d
  Delta <- as.list(V)
  for (k in 2:d){
    toremove <- list()
    Deltai <- list()
    for (C in utils::combn(V, k, simplify = FALSE)){
      if ((utils::combn(C, k-1, simplify = FALSE) %in% Delta) && is.M(solve(S[C,C]))) {
        Deltai <- append(Deltai,list(C))
        toremove <- append(toremove,utils::combn(C, k-1, simplify = FALSE))}
    }
    Delta <- setdiff(Delta,toremove)
    Delta <- append(Delta,Deltai)
    if (length(Deltai)<=k+1){
      return(Delta)
      break
    }
  }
}
