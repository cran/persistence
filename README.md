
<!-- README.md is generated from README.Rmd. Please edit that file -->

# persistence

<!-- badges: start -->

<!-- badges: end -->

> **Deprecated.** `persistence` has been superseded by the
> [`scalednap`](https://CRAN.R-project.org/package=scalednap) package — a strict
> superset with the same functions and more — and is being retired from CRAN.
> Please install `scalednap` instead: `install.packages("scalednap")` (the same
> algorithm is also available for Python: `pip install scalednap`).

The goal of persistence is to …

## Installation

You can install the stable version of **persistence** from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("persistence")
```

If you want to try the latest development version with the newest
updates, you can install it from [CRAN](https://CRAN.R-project.org)
using:

``` r
# install.packages("pak")
pak::pak("aavellone/persistenceR")
```

## Example

This is a basic example which shows you how to solve a common problem
using **persistence**:

``` r
library(persistence)

# --- EXAMPLE 1: Standard input (vectors and matrices) ---
edg <- c(1, 2, 1, 3, 1, 4, 2, 3, 3, 4, 4, 5, 5, 6, 5, 7, 6, 7)
edge_list <- matrix(edg, ncol = 2, byrow = TRUE)
vertex <- c(1, 2, 3, 4, 5, 6, 7)

cluster_milano(x = vertex, edge_list = edge_list)
#> $membership
#> [1] 1 1 1 1 2 2 2
#> 
#> $score
#> [1] 0.7662338
#> 
#> $seed
#> [1] "3563829245"

# --- EXAMPLE 2: igraph input ---
if (requireNamespace("igraph", quietly = TRUE)) {
  g <- igraph::make_ring(10)
  cluster_milano(g)
}
#> $membership
#>  [1] 1 1 2 2 3 3 4 4 5 5
#> 
#> $score
#> [1] 1.5
#> 
#> $seed
#> [1] "806094222"
```
