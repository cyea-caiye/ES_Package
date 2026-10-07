#' Printing Summaries of Autoregressive Quantile-Expected Shortfall Fits
#'
#' @description
#' \code{print.summary.estools.ts} prints a summary of inference results from
#' a fitted autoregressive quantile-expected shortfall model.
#'
#' @param x an object of class \code{"summary.estools.ts"}.
#' @param digits number of significant digits used for printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The printed output includes the matched function call, quantile level,
#' conditional quantile and expected shortfall model specifications,
#' confidence level, and the inference table for the final weighted expected
#' shortfall coefficients.
#'
#' @return
#' Invisibly returns \code{x}.
#'
#' @seealso
#' \code{\link{summary.estools.ts}},
#' \code{\link{estools.ts}}
#'
#' @export
print.summary.estools.ts <- function(
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
