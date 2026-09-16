#' @rdname summary.lesreg
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
