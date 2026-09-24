#' @rdname lesreg.highdim.inference
#' @param digits number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#' @method print lesreg.highdim.inference
#' @export
print.lesreg.highdim.inference <- function(x, digits = max(5L, getOption("digits") - 2L), ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
  }
  cat("\nHigh-dimensional expected shortfall inference\n")
  if(!is.null(x$col)){
    cat("\nCoefficient index:", x$col, "\n")
  }
  if(!is.null(x$tau)){
    cat("Tau:", x$tau, "\n")
  }
  if(!is.null(x$tail)){
    cat("Tail:", x$tail, "\n")
  }
  if(!is.null(x$variance.method)){
    cat("Variance method:", x$variance.method, "\n")
  }
  if(!is.null(x$conf.level)){
    cat("Confidence level:", format(round(x$conf.level, digits = digits)), "\n")
  }
  cat("\nDebiased estimate:\n")
  print(format(round(x$theta_debias, digits = digits)), quote = FALSE)
  cat("\nConfidence interval:\n")
  print(format(round(x$Conf.int, digits = digits)), quote = FALSE)
  invisible(x)
}
