#' Computes the Z-matrix
#'
#' This function oututs the Z matrix, that is, the unique dominating ultrametric.
#' @param S a covariance matrix
#' @keywords xxx
#' @export
#' @examples
#' 
Zmatrix <- function(S){
  p <- nrow(S)
  R <- cov2cor(S)
  # Compute the distances. Non-positive correlations correspond to very big distances.
  D <- as.dist(-log((R>0)*R+(R<=0)*1e-20))
  # use single-linkage clustering method in R
  hcl <- hclust(D,method="single")
  # recover how hclust merges variables and use it to recover the corresponding ultrametric 
  subs <- list()
  length(subs) <- p-1
  Z <- matrix(0,p,p)
    for (i in 1:(p-1)){
      if (hcl$merge[i,1]<0 && hcl$merge[i,2]<0) {
        subs[[i]] <- union(-hcl$merge[i,1],-hcl$merge[i,2])
        Z[-hcl$merge[i,1],-hcl$merge[i,2]]<- Z[-hcl$merge[i,2],-hcl$merge[i,1]] <- hcl$height[i]}
      else if (hcl$merge[i,1]<0 && hcl$merge[i,2]>0) {
        subs[[i]] <- union(-hcl$merge[i,1],subs[[hcl$merge[i,2]]])
        Z[-hcl$merge[i,1],subs[[hcl$merge[i,2]]]] <- hcl$height[i]
        Z[subs[[hcl$merge[i,2]]],-hcl$merge[i,1]] <- hcl$height[i]
      }
      else {
        subs[[i]] <- union(subs[[hcl$merge[i,1]]],subs[[hcl$merge[i,2]]])
        Z[subs[[hcl$merge[i,1]]],subs[[hcl$merge[i,2]]]] <- hcl$height[i]
        Z[subs[[hcl$merge[i,2]]],subs[[hcl$merge[i,1]]]] <- hcl$height[i]
          }
    }
  # Z is the corresponding ultrametric. Now compute the indiced correlation matrix
  Z <- exp(-Z)
  # finally output the corresponding covariance matrix
  return(diag(sqrt(diag(S)))%*%Z%*%diag(sqrt(diag(S))))
}

