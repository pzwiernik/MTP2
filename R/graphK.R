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
graphK <- function(S,K,names=FALSE,tol=1e-10,root=c()){
  p <- nrow(S)
  if (names==FALSE){
    dimnames(S) <- dimnames(K) <- NULL
  }
  # the MLE graph
  G <- igraph::graph_from_adjacency_matrix(1*(K< -tol),mode="undirected")
  # the MWSF
  T <- MWSF(S,positive=TRUE)
  V(T)$name <- V(G)$name
  # the excessive correlation graph
  C <- CLmatrix(S,positive=TRUE)
  edges <- c()
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      if (S[i,j]-C[i,j]>= -tol) edges <- c(edges,i,j)
    }
  }
  ecg <- igraph::make_empty_graph(p) 
  ecg <- igraph::as.undirected(igraph::add_edges(ecg,edges))
  V(ecg)$name <- V(G)$name
  
  # check if the graphs are nested as predicted by the theory (numerical issues still possible)
  if (igraph::ecount(igraph::difference(T,G))>0 || igraph::ecount(igraph::difference(G,ecg))>0) {
    print("Something is wrong. The three graphs are not nested.")
  }

  
  # format the output 
  edg <- igraph::get.edges(ecg,igraph::E(ecg))
  Gedges <- igraph::get.edges(G,igraph::E(G))
  Tedges <- igraph::get.edges(T,igraph::E(T))
  igraph::E(ecg)$color <- "#C9C9C9"
  igraph::E(ecg)$width<- 1
  for (i in 1:nrow(edg)){
    if (TRUE %in% apply(Tedges, 1, function(x, want) isTRUE(all.equal(x, want)), edg[i,])) {
      igraph::E(ecg)[i]$color <- "#D53E4F"
      igraph::E(ecg)[i]$width <- 4
    }
    else if (TRUE %in% apply(Gedges, 1, function(x, want) isTRUE(all.equal(x, want)), edg[i,])){
      igraph::E(ecg)[i]$color <- "#3288BD"
      igraph::E(ecg)[i]$width <- 2
    }
  }
  igraph::V(ecg)$size <- 7
  

  #plot the graph
  if (sum(root)==0) {
    igraph::plot.igraph(ecg,layout=igraph::layout_as_tree(T,circular = TRUE))
  }
  if (sum(root)>0) {
    igraph::plot.igraph(ecg,vertex.label=1:p,layout=igraph::layout_as_tree(T,root=root,circular = TRUE))
  }
  #igraph::tkplot(ecg)
  }
