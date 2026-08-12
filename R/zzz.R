#' @importFrom methods hasArg
#' @importFrom utils packageVersion combn head
#' @importFrom stats sd var
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @import RcppParallel

NULL

# Null-coalescing operator for R < 4.4 compatibility
if(!exists("%||%", baseenv())) {
  `%||%` = function(x, y) if(is.null(x)) y else x
}
