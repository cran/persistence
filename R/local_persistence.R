#' @title local_persistence
#'
#' @description
#' Given the incidence vector of a vertex subset, calculates either the
#' persistence probability or the null-adjusted persistence of cluster C.
#' This function is polymorphic: it automatically detects the input type and accepts
#' either a vertex vector (accompanied by an edge list) or directly an \code{igraph} object.
#'
#' @param x   An integer or character vector representing the graph vertices,
#'            OR an object of class \code{igraph}.
#' @param ... Additional arguments passed to specific methods (e.g., \code{edge_list},
#'            \code{weights}, \code{cluster}, \code{H0}).
#'
#' @return Numeric scalar: the persistence probability when \code{H0 = NULL},
#'         the null-adjusted persistence when \code{H0 = 0},
#'         or the null-adjusted persistence density when \code{H0} is in \eqn{(0, 1]}.
#'
#' @examples
#' library(persistence)
#'
#' # --- EXAMPLE 1: Standard input (vectors and matrices) ---
#' edg <- c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 6, 7)
#' edge_list <- matrix(edg, ncol = 2, byrow = TRUE)
#' vertex <- c(1, 2, 3, 4, 5, 6, 7)
#' cluster_bin <- c(1, 1, 1, 1, 0, 0, 0)
#' local_persistence(x = vertex, edge_list = edge_list, cluster = cluster_bin)
#'
#' # --- EXAMPLE 2: igraph input ---
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   g <- igraph::make_ring(10)
#'   cluster_bin <- c(rep(1, 5), rep(0, 5))
#'   local_persistence(g, cluster = cluster_bin)
#' }
#'
#' @export
local_persistence <- function(x, ...) {
  UseMethod("local_persistence")
}

#' @rdname local_persistence
#' @param edge_list  Integer matrix with two columns representing the graph edge list.
#' @param weights    Numeric vector of edge weights. If \code{NULL}, all weights default to 1.
#' @param cluster    Binary incidence vector of the cluster: \code{x_i = 1} if \code{i in C}, 0 otherwise.
#' @param H0         Optional numeric value in \eqn{[0, 1]}, or \code{NULL}. Default is \code{0}.
#' @export
local_persistence.default <- function(x, edge_list, weights = NULL, cluster, H0 = 0, ...) {
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

  # --- Validate Cluster -------------------------------------------------------
  if (is.logical(cluster)) {
    cluster <- as.integer(cluster)
  }
  if (!is.numeric(cluster)) {
    stop("cluster must be a binary (0/1) vector.")
  }
  if (length(cluster) != length(vertex)) {
    stop(sprintf("cluster must have length %d.", length(vertex)))
  }
  if (!all(cluster %in% c(0, 1))) {
    stop("cluster must contain only 0 and 1.")
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

  if (!is.integer(cluster)) {
    cluster <- as.integer(cluster)
  }

  result <- tryCatch(
    .Call(C_local_persistence, vertex, edge_list, weights, cluster, H0),
    error = function(e) stop(e$message)
  )

  return(result)
}

#' @rdname local_persistence
#' @export
local_persistence.igraph <- function(x, cluster, H0 = 0, ...) {
  g_data <- .extract_from_igraph(x)

  local_persistence.default(
    x = g_data$vertex,
    edge_list = g_data$edge_list,
    weights = g_data$weights,
    cluster = cluster,
    H0 = H0
  )
}
