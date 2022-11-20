#' evaluate additional expression inide pipe syntax.
#'
#' @param x placeolder for actual data passing througt the pipe
#' @param expr expression to evaluate at step usually message function
#' @return result is unaltered input from x

log_midpipe <- function(x, expr) {
  result <- x
  expr 
  result
}