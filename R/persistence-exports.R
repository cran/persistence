#' @title Persistence (deprecated)
#'
#' @description
#' Given a non-oriented graph, calculates the optimal vertex partition using
#' persistence as the objective function.
#'
#' \strong{This package is deprecated} and has been superseded by the
#' \pkg{scalednap} package, a strict superset that provides the same functions
#' and additional functionality. Please install \pkg{scalednap} instead
#' (\code{install.packages("scalednap")}); the same algorithm is also available
#' for Python (\code{pip install scalednap}).
#'
#' @details
#' See manual entries.
#'
#' @docType package
#' @useDynLib persistence, .registration=TRUE, .fixes="C_"
"_PACKAGE"
