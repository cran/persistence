.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "'persistence' is deprecated and is being retired from CRAN.\n",
    "It has been superseded by the 'scalednap' package, a strict superset ",
    "with the same functions and more.\n",
    "Please switch to 'scalednap':\n",
    "    install.packages(\"scalednap\")\n",
    "The same algorithm is also available for Python: pip install scalednap"
  )
}
