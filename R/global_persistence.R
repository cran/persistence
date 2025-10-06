#' @title
#' global persistence
#'
#' @description
#' Given a partition of the graph vertices, it calculates the global persistence as the sum of the persistences
#' of the single clusters. Persistence can be referred to the null-adjusted o to the probability.
#'
#' @param vertex
#' the vertices of the graph, whose label are integers and they must be consistent with the edge sets.
#'
#' @param edge_list
#' the graph edge list in the form of an integer matrix with two columns.
#'
#' @param weights
#' the graph edge weights.
#'
#' @param membership
#' An integer vector representing the vertex membership:  x_i = k if i in C_k.
#'
#' @param H0
#' If true, it calculates the null-adjusted persistence, if false, the persistence probability.
#'
#' @returns value
#' A list containing the following:
#' \describe{
#'  \item{value}{The global persistence of the partition.}
#'  \item{clusters_value}{The local persistence of each cluster. If for some k we have v_k = NaN, then
#'  C_k is empty in the input membership.}
#' }
#'
#'
#' @examples
#' library(persistence)
#' library(igraph)
#'
#' edg = c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 7, 9, 8, 9)
#' vertex = unique(edg)
#' edg = t(matrix(as.integer(edg), nrow = 2 ))
#' rete <- graph_from_edgelist(edg, directed = FALSE) # I graph this matrix
#' plot(rete)
#'
#' membership = c(1, 1, 1, 1, 2, 2, 2, 2, 2)
#' v1 = global_persistence(vertex, edg, weights=NULL, membership, H0=TRUE)
#' print(paste("global null-adjusted persistence: ", v1$value))
#' print(paste("null-adjusted persistence per cluster: ", v1$clusters_value))

#' @name global_persistence
#' @export global_persistence
global_persistence <- function(vertex, edge_list, weights, membership, H0 = TRUE) {
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

  if (!is.numeric(membership)) {
    stop("membership must be an array of integer")
  }

  if (length(membership) != length(vertex)) {
    stop(paste("membership must be of length:", length(vertex), sep=""))
  }

  if (!is.logical(H0) || length(H0) != 1) {
    stop("H0 must be TRUE/FALSE")
  }

  tryCatch({
    result <- .Call("_global_persistence", vertex, edge_list, weights, as.integer(membership), as.logical(H0))
  }, warning = function(war) {
    # warning handler picks up where error was generated
    print(paste("MY_WARNING:  ",war))
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("MY_ERROR:  ",err))
  }) # END tryCatch
  return (result)
}




