#' @title
#' local_persistence
#'
#' @description
#' Given the incidence vector of a vertex subset, it calculates the persistence probability
#' or the null-adjusted persistence of C.
#' @param vertex
#' the vertices of the graph, whose label are integers and they must be consistent with the edge sets
#'
#' @param edge_list
#' the graph edge list in the form of an integer matrix with two columns
#'
#' @param weights
#' the graph edge weights.
#'
#' @param cluster
#' A binary vector representing the incidence vector of the cluster:  x_i = 1 if i in C, 0 otherwise.
#'
#' @param H0
#' if true, it calculates the null-adjusted persistence, if false, the persistence probability.
#'
#' @return
#' the value of the null-adjusted persistence if H0 = T, the value of the persistence probability if H0 = F
#'
#' @examples
#' library(persistence)
#' library(igraph)
#'
#' edg = c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 7, 9, 8, 9)
#' print(length(edg) / 2.0)
#' vertex = unique(edg)
#' edg = t(matrix(as.integer(edg), nrow = 2 ))
#' rete <- graph_from_edgelist(edg, directed = FALSE) # I graph this matrix
#' plot(rete)
#'
#' cluster = rep(0, length(vertex))
#' v1 = c(1, 2, 3, 4)
#' cluster[v1] = 1
#' f1 = local_persistence(vertex, edg, weights=NULL, cluster, H0 = TRUE)
#' f2 = local_persistence(vertex, edg, weights=NULL, cluster, H0 = FALSE)
#'
#' @name local_persistence
#' @export local_persistence
local_persistence <- function(vertex, edge_list, weights, cluster, H0 = TRUE) { #vettore 0/1
  if (!is.vector(vertex)) {
    stop("vertex must be an array")
  }
  if (length(unique(vertex)) != length(vertex)) {
    stop("vertex contains duplicated values")
  }
  if (!is.matrix(edge_list)) {
    stop("edge_list must be an edge list 1")
  }

  if (ncol(edge_list) != 2) {
    stop("edge_list must be an edge list 2")
  }
  if (!(all(unique(as.vector(edge_list)) %in% vertex))) {
    stop("edge_list contains values not belonging to vertex")
  }

  vertex <- as.character(vertex)
  edge_list <- matrix(as.character(edge_list), ncol = 2)

  if (is.null(weights)) {
    weights <- as.numeric(rep(1.0, nrow(edge_list)))
  }

  if (!is.vector(weights)) {
    stop("weights must be an array")
  }

  if (length(weights) != nrow(edge_list)) {
    stop(paste("weights must be of length:", nrow(edge_list), sep=""))
  }

  if (is.logical(cluster)) {
    cluster <- ifelse(cluster == TRUE, 1, 0)
  }
  if (!is.numeric(cluster)) {
    stop("cluster must be an array of 0/1")
  }
  if (length(unique(as.vector((edge_list)))) != length(cluster)) {
    stop(paste("cluster must be of length:", length(vertex), sep=""))
  }
  if (length(setdiff(unique(cluster), c(1, 0))) != 0) {
    stop("cluster must be an array of 0/1")
  }
  if (!is.logical(H0) || length(H0) != 1) {
    stop("H0 must be TRUE/FALSE")
  }
  tryCatch({
    result <- .Call("_local_persistence", vertex, edge_list, weights, as.integer(cluster), as.logical(H0))
  }, warning = function(war) {
    # warning handler picks up where error was generated
    print(paste("MY_WARNING:  ",war))
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("MY_ERROR:  ",err))
  }) # END tryCatch
  return (result)
}
