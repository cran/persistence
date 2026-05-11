#' @title Extract components from an igraph object
#' @description Internal utility function to extract vertices, edges, and weights.
#' @noRd
.extract_from_igraph <- function(x) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("The 'igraph' package is required for this operation. Install it with install.packages('igraph').")
  }
  if (!igraph::is_igraph(x)) {
    stop("'x' must be a valid igraph object.")
  }
  if (igraph::is_directed(x)) {
    stop("The graph 'x' must be undirected.")
  }

  weights <- igraph::E(x)$weight

  if (is.null(igraph::V(x)$name)) {
    vertex <- as.character(as.numeric(igraph::V(x)))
    edge_list <- igraph::as_edgelist(x, names = FALSE)
    mode(edge_list) <- "character"
  } else {
    vertex <- igraph::V(x)$name
    edge_list <- igraph::as_edgelist(x, names = TRUE)
  }

  list(vertex = vertex, edge_list = edge_list, weights = weights)
}

#' @title cluster_milano
#'
#' @description
#' Calculates the vertex partition with maximum global null-adjusted persistence.
#' This function is polymorphic: it automatically detects the input type and accepts
#' either a vertex vector (accompanied by an edge list) or directly an \code{igraph} object.
#'
#' @param x   An integer or character vector representing the graph vertices,
#'            OR an object of class \code{igraph}.
#' @param ... Additional arguments passed to specific methods (e.g., \code{edge_list},
#'            \code{weights}, \code{tol}, \code{max_level}, etc.).
#'
#' @return A list with three elements:
#' \describe{
#'   \item{membership}{The optimal vertex partition.}
#'   \item{score}{The measure value of the optimal partition.}
#'   \item{seed}{The seed used to generate random numbers.}
#' }
#'
#' @examples
#' library(persistence)
#'
#' # --- EXAMPLE 1: Standard input (vectors and matrices) ---
#' edg <- c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 6, 7)
#' edge_list <- matrix(edg, ncol = 2, byrow = TRUE)
#' vertex <- c(1, 2, 3, 4, 5, 6, 7)
#' cluster_milano(x = vertex, edge_list = edge_list)
#'
#' # --- EXAMPLE 2: igraph input ---
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   g <- igraph::make_ring(10)
#'   cluster_milano(g)
#' }
#'
#' @export
cluster_milano <- function(x, ...) {
  UseMethod("cluster_milano")
}

#' @rdname cluster_milano
#' @param edge_list  Integer matrix with two columns representing the graph edge list.
#' @param weights    Numeric vector of positive edge weights. If \code{NULL}, all weights default to 1.
#' @param membership Integer vector representing the starting partition: \code{x_i = k} if \code{i in C_k}.
#'                   If \code{NULL}, each vertex starts in its own cluster.
#' @param H0         Logical value. Default is \code{TRUE}.
#'                   If \code{TRUE}, returns the null-adjusted persistence.
#'                   If \code{FALSE}, returns the persistence probability.
#' @param seed       Non-negative integer seed for the random number generator.
#'                   If \code{NULL}, an internal default is used.
#' @param tol        Optional numeric tolerance for the stopping criterion.
#'                   If \code{NULL} (default), an adaptive threshold is calculated dynamically in C++.
#' @param max_level  Optional integer representing the maximum number of aggregation levels.
#'                   If \code{0} (default) or \code{NULL}, the algorithm runs until convergence.
#' @export
cluster_milano.default <- function(x, edge_list, weights = NULL, membership = NULL, H0 = TRUE, seed = NULL, tol = NULL, max_level = 0L, ...) {
  vertex <- x
  n <- length(vertex)

  # --- Validate vertex and edge_list ------------------------------------------
  if (!is.vector(vertex)) {
    stop("vertex must be a vector.")
  }
  if (any(duplicated(vertex))) {
    stop("vertex contains duplicated values.")
  }
  if (!is.matrix(edge_list) || ncol(edge_list) != 2L) {
    stop("edge_list must be a two-column integer matrix.")
  }
  if (!all(unique(as.vector(edge_list)) %in% vertex)) {
    stop("edge_list contains values not belonging to vertex.")
  }
  # --- Validate weights -------------------------------------------------------
  if (!is.null(weights)) {
    if (!is.vector(weights)) stop("weights must be a vector.")
    if (length(weights) != nrow(edge_list)) stop(sprintf("weights must have length %d.", nrow(edge_list)))
  }
  # --- Validate Membership ----------------------------------------------------
  if (is.null(membership)) {
    membership <- seq_len(n)
  } else {
    if (length(membership) != n) {
      stop("membership and vertex must agree on length.")
    }
    if (min(membership) < 1L || max(membership) > n) {
      stop(sprintf("membership values must be between 1 and %d.", n))
    }
  }
  # --- Validate H0 ------------------------------------------------------------
  if (!is.logical(H0) || length(H0) != 1L || is.na(H0)) {
    stop("H0 must be a single logical value (TRUE or FALSE).")
  }
  H0_internal <- if (H0) 0 else NULL
  # --- Validate Seed ----------------------------------------------------------
  if (!is.null(seed)) {
    if (length(seed) != 1L || seed < 0 || seed != round(seed)) {
      stop("seed must be a single non-negative integer.")
    }
  }
  # --- Validate tol -----------------------------------------------------------
  if (!is.null(tol)) {
    if (!is.numeric(tol) || length(tol) != 1L || tol <= 0) {
      stop("tol must be a single positive numeric value, or NULL.")
    }
  }
  # --- Validate max_level -----------------------------------------------------
  if (is.null(max_level)) {
    max_level <- 0L
  }
  if (!is.numeric(max_level) || length(max_level) != 1L || max_level < 0 || max_level != round(max_level)) {
    stop("max_level must be a single non-negative integer.")
  }
  if (!is.character(vertex)) vertex <- as.character(vertex)
  if (!is.character(edge_list)) edge_list <- matrix(as.character(edge_list), ncol = 2L)

  if (is.null(weights)) {
    weights <- rep(1.0, nrow(edge_list))
  } else if (!is.numeric(weights)) {
    weights <- as.numeric(weights)
  }

  if (is.null(membership)) {
    membership <- seq_len(n)
  } else if (!is.integer(membership)) {
    membership <- as.integer(membership)
  }

  if (!is.null(seed) && !is.integer(seed)) seed <- as.integer(seed)
  if (!is.null(tol) && !is.double(tol)) tol <- as.numeric(tol)
  if (!is.integer(max_level)) max_level <- as.integer(max_level)

  result <- tryCatch(
    .Call(C_cluster_milano, vertex, edge_list, weights, membership, H0_internal, seed, tol, max_level),
    error = function(e) stop(e$message)
  )

  return(result)
}

#' @rdname cluster_milano
#' @export
cluster_milano.igraph <- function(x, membership = NULL, H0 = TRUE, seed = NULL, tol = NULL, max_level = 0L, ...) {
  g_data <- .extract_from_igraph(x)

  cluster_milano.default(
    x = g_data$vertex,
    edge_list = g_data$edge_list,
    weights = g_data$weights,
    membership = membership,
    H0 = H0,
    seed = seed,
    tol = tol,
    max_level = max_level
  )
}
