#' @rdname summary.lesreg.ts
#' @param x an object of class \code{"summary.lesreg.ts"}.
#' @param digits number of significant digits used for printing.
#' @export
print.summary.lesreg.ts <- function(
    x,
    digits = max(5L, getOption("digits") - 2L),
    ...){
  cat("\nCall: ")
  dput(x$call)
  cat("\nTau: ")
  print(format(round(x$tau, digits = digits)), quote = FALSE)
  cat("Quantile model: ", x$qmodel, "\n", sep = "")
  cat("Expected shortfall model: ", x$emodel, "\n", sep = "")
  cat("Confidence level: ")
  print(format(round(x$level, digits = digits)), quote = FALSE)
  cat("\nWeighted expected shortfall coefficients:\n")
  print(format(round(x$coefficients, digits = digits)), quote = FALSE)
  invisible(x)
}
