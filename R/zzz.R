#' @title Start up message
#' @importFrom utils citation
#' @keywords internal
#' @noRd
.onAttach <- function(...) {
  packageStartupMessage("\nType 'citation(\"TestDimorph\")' for citing this R package in publications.")
}
