#' @title
#' cluster Milano
#'
#' @description
#' Calculates the partition with maximum global null-adjuted persistence.
#'
#' @param vertex
#' the vertices of the graph, whose label are integers and they must be consistent with the edge sets.
#'
#'
#' @param edge_list
#' the graph edge list in the form of an integer matrix with two columns.
#'
#' @param seed
#' As some steps of the algorithm are random, users may experiments with different seeds of random numbers.
#' @return
#' A list containing:
#' \describe{
#' \item{membeship}{The optimal vertex partition.}
#' \item{value}{The null-adjusted persistence of the partition.}
#' \item{seed}{The used seed to generate random numbers.}
#' }
#'
#' @examples
#' library(persistence)
#' library(igraph)
#'
#' edg = c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 7, 9, 8, 9)
#' print(length(edg) / 2.0)
#' vertex = unique(edg)
#' edg = t(matrix(as.integer(edg), nrow = 2 ))
#' rete <- graph_from_edgelist(edg, directed = FALSE)
#' plot(rete)
#' seed <- sample(1:as.integer(.Machine$integer.max),1, replace= FALSE)
#' r = cluster_milano(vertex, edg, seed=seed)
#' print(paste("The optimal null-adjusted persistence is: ", r$measure))
#' print(paste("The optimal persistence probability is: ", r$measure + 1))
#'
#'
#' @name cluster_milano
#' @export cluster_milano
cluster_milano <- function(vertex, edge_list, seed=NULL) {
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

  if (!is.null(seed) && (seed < 0 || seed != round(seed))) {
    stop("seed must be a positive integer")
  }
  seed <- as.integer(seed)
  result <- NULL
  tryCatch({
    result <- .Call("_cluster_milano", vertex, edge_list, seed)
  }, warning = function(war) {
    # warning handler picks up where error was generated
    print(paste("MY_WARNING:  ",war))
  }, error = function(err) {
    # error handler picks up where error was generated
    print(paste("MY_ERROR:  ",err))
  }) # END tryCatch
  return (result)
}
