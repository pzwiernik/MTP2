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
viol3 <- function(S){
    p <- nrow(S)
    vio <- rep(0,p)
    for (i in 1:(p-2)){
      for (j in (i+1):(p-1)){
        for (k in (j+1):p){ 
          vio[i] <- vio[i]+min(c(S[i,j]*S[i,k]*S[j,k],0))^2
          vio[j] <- vio[j]+min(c(S[i,j]*S[i,k]*S[j,k],0))^2
          vio[k] <- vio[k]+min(c(S[i,j]*S[i,k]*S[j,k],0))^2
        } 
      }
    }
    return(vio)
}
  
