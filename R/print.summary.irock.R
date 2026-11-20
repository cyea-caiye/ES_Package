#' @rdname summary.irock
#' @export
print.summary.irock <- function(x, digits = max(5L, getOption("digits") - 2L),
                                ...) {
  if (!is.null(x$call)) {
    cat("\nCall: ")
    dput(x$call)
  }
  cat("\nTail: ", x$tail, "\n", sep = "")
  cat("Tau: ", format(x$tau, digits = digits), "\n", sep = "")
  cat("Variance method: ", x$method, " (B = ", x$B, ")\n", sep = "")
  cat("Confidence level: ", format(x$level, digits = digits), "\n", sep = "")
  coefficients <- x$coefficients.results
  p <- coefficients[, "P-value"]
  cutoff <- 10^(-digits)
  small <- !is.na(p) & p < cutoff
  coefficients[small, "P-value"] <- NA_real_
  display <- coefficients
  for (j in seq_len(ncol(coefficients))) {
    display[, j] <- format(coefficients[, j], digits = digits,
                           scientific = FALSE, trim = TRUE)
  }
  display[small, "P-value"] <- paste0("< ", format(cutoff, scientific = FALSE))
  cat("\nExpected shortfall regression coefficients:\n")
  print(display, quote = FALSE, ...)
  invisible(x)
}
