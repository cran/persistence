#' @title global_persistence
#'
#' @description
#' Given a partition of the graph vertices, calculates the global persistence
#' as the sum of the local persistences of the individual clusters.
#' Persistence can be either null-adjusted or probability-based.
#' This function is polymorphic: it automatically detects the input type and accepts
#' either a vertex vector (accompanied by an edge list) or directly an \code{igraph} object.
#'
#' @param x   An integer or character vector representing the graph vertices,
#'            OR an object of class \code{igraph}.
#' @param ... Additional arguments passed to specific methods (e.g., \code{edge_list},
#'            \code{weights}, \code{membership}, \code{H0}).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{score}{The global persistence of the partition.}
#'   \item{clusters_value}{The local persistence of each cluster. A value of \code{NaN}
#'     indicates that cluster \eqn{C_k} is empty in the input \code{membership}.}
#' }
#'
#' @examples
#' library(persistence)
#'
#' # --- EXAMPLE 1: Standard input (vectors and matrices) ---
#' edg <- c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 6, 7)
#' edge_list <- matrix(edg, ncol = 2, byrow = TRUE)
#' vertex <- c(1, 2, 3, 4, 5, 6, 7)
#' mem <- c(1, 1, 1, 1, 2, 2, 2)
#' global_persistence(x = vertex, edge_list = edge_list, membership = mem)
#'
#' # --- EXAMPLE 2: igraph input ---
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   g <- igraph::make_ring(10)
#'   mem <- c(rep(1, 5), rep(2, 5))
#'   global_persistence(g, membership = mem)
#' }
#'
#' @export
global_persistence <- function(x, ...) {
  UseMethod("global_persistence")
}

#' @rdname global_persistence
#' @param edge_list  Integer matrix with two columns representing the graph edge list.
#' @param weights    Numeric vector of edge weights. If \code{NULL}, all weights default to 1.
#' @param membership Integer vector of vertex cluster assignments: \code{x_i = k} if \code{i in C_k}.
#' @param H0         Optional numeric value in \eqn{[0, 1]}, or \code{NULL}. Default is \code{0}.
#' @export
global_persistence.default <- function(x, edge_list, weights = NULL, membership, H0 = 0, ...) {
  vertex <- x

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

  # --- Validate H0 ------------------------------------------------------------
  if (!is.null(H0)) {
    if (!is.numeric(H0) || length(H0) != 1L) {
      stop("H0 must be a single numeric value between 0 and 1, or NULL.")
    }
    if (H0 < 0 || H0 > 1) {
      stop("H0 must be between 0 and 1.")
    }
  }

  # --- Coercions (Zero-Copy) --------------------------------------------------
  if (!is.character(vertex)) {
    vertex <- as.character(vertex)
  }
  if (!is.character(edge_list)) {
    edge_list <- matrix(as.character(edge_list), ncol = 2L)
  }

  if (is.null(weights)) {
    weights <- rep(1.0, nrow(edge_list))
  } else {
    if (!is.vector(weights)) {
      stop("weights must be a vector.")
    }
    if (length(weights) != nrow(edge_list)) {
      stop(sprintf("weights must have length %d.", nrow(edge_list)))
    }
    if (!is.numeric(weights)) {
      weights <- as.numeric(weights)
    }
  }

  if (length(membership) != length(vertex)) {
    stop("membership and vertex must agree on length.")
  }
  if (!is.integer(membership)) {
    membership <- as.integer(membership)
  }

  # --- C call (RIPRISTINATA) --------------------------------------------------
  result <- tryCatch(
    .Call(C_global_persistence, vertex, edge_list, weights, membership, H0),
    error = function(e) stop(e$message)
  )

  return(result)
}

#' @rdname global_persistence
#' @export
global_persistence.igraph <- function(x, membership, H0 = 0, ...) {
  g_data <- .extract_from_igraph(x)

  global_persistence.default(
    x = g_data$vertex,
    edge_list = g_data$edge_list,
    weights = g_data$weights,
    membership = membership,
    H0 = H0
  )
}
