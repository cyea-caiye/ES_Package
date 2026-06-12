#' Print Linear Expected Shortfall Regression Summary Object
#'
#' @description
#' Print summary of a linear expected shortfall regression object.
#'
#' @usage
#' \method{print}{summary.lesreg}(x, digits = max(5L, getOption("digits") - 2L), ...)
#'
#' @param x An object of class \code{"summary.lesreg"}, typically produced by
#'   a call to \code{\link{summary.lesreg}}.
#' @param digits Significant digits reported in the printed table.
#' @param ... Optional arguments passed to the printing function.
#'
#' @seealso
#' \code{\link{summary.lesreg}}, \code{\link{lesreg}}
#'
#' @export

print.summary.lesreg <- function(x, digits = max(5L, getOption("digits") - 2L), ...){
  cat("\nCall: ")
  dput(x$call)
  cat("\nTail: ")
  cat(x$tail, "\n")
  cat("\nTau: ")
  print(format(round(x$tau, digits = digits)), quote = FALSE, ...)
  cat("Variance method: ")
  cat(x$method, "\n")
  cat("Confidence level: ")
  print(format(round(x$level, digits = digits)), quote = FALSE, ...)
  cat("\nCoefficients:\n")
  print(format(round(x$coefficients, digits = digits)), quote = FALSE)
  invisible(x)
}
