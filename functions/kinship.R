
find_maximum_independet_set <- function(kinships) {
  
  oldval <- igraph::igraph_opt("add.vertex.names")
  on.exit(igraph::igraph_options(add.vertex.names = oldval), add = TRUE)
  igraph::igraph_options(add.vertex.names = FALSE)
  
  graph <- igraph::graph_from_data_frame(d = kinships[, .(ID1, ID2)], directed = FALSE)
  
  S <- character(0)
  
  total_length <- length(igraph::V(graph))
  
  pb <- progress::progress_bar$new(
    format = " [:bar] :current/:total (:percent) :eta",
    total = total_length
  )
  
  pb$tick(0)
  
  while (l <- length(vg <- igraph::V(graph))) {
    mv <- which.min(.Call(igraph:::C_R_igraph_degree, graph, 1:l-1, 1, FALSE)) # find vertex with minimum degree, i.e. new independent vertex
    .Call(igraph:::C_R_igraph_finalizer)
    S <- c(S, vg$name[mv]) # add independent vertex name to set of independent samples
    R <- c(mv, .Call(igraph:::C_R_igraph_neighbors, graph, mv-1, 3)+1) # find all neighbors of new independent vertex
    .Call(igraph:::C_R_igraph_finalizer)
    graph <- .Call(igraph:::C_R_igraph_delete_vertices, graph, R-1) # remove new independent vertex and all its neighbors from graph
    .Call(igraph:::C_R_igraph_finalizer)
    pb$tick(length(R))
  }
  
  return(S)
  
}
