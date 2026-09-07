# persistence 1.0.1

* **Deprecation.** `persistence` has been superseded by the
  [`scalednap`](https://CRAN.R-project.org/package=scalednap) package — a strict
  superset that provides the same functions (`cluster_milano()`,
  `global_persistence()`, `local_persistence()`) and additional functionality —
  and is being retired from CRAN. Please switch to `scalednap`
  (`install.packages("scalednap")`; the same algorithm is also available for
  Python via `pip install scalednap`). A startup message now points users to the
  replacement. This is the final release of `persistence`.
