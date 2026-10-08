#' @rdname summary.esreg.ts
#' @export
print.summary.esreg.ts <- function(
    x,
    digits = max(5L, getOption("digits") - 2L),
    ...){
  cat("\nCall: ")
  dput(x$call)
  cat("\nTau: ")
  print(format(round(x$tau, digits = digits)), quote = FALSE)
  cat("Quantile model: ", x$qmodel, "\n", sep = "")
  cat("Expected shortfall model: ", x$emodel, "\n", sep = "")
  cat("Covariance method: ", x$covMethod, "\n", sep = "")
  cat("Truncation method: ", x$truncationMethod, "\n", sep = "")
  cat("Confidence level: ")
  print(format(round(x$level, digits = digits)), quote = FALSE)
  cat("\nWeighted expected shortfall coefficients:\n")
  print(format(round(x$coefficients, digits = digits)), quote = FALSE)
  invisible(x)
}
