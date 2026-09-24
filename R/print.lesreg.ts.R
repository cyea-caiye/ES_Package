#' Print an Autoregressive Quantile-Expected Shortfall Fit
#'
#' @param x An object of class \code{"lesreg.ts"}.
#' @param digits Number of significant digits used for printing.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.lesreg.ts <- function(
    x,
    digits = max(5L, getOption("digits") - 2L),
    ...){
  if(!inherits(x, "lesreg.ts")){
    stop("'x' must be a lesreg.ts object")
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\nTau:", x$tau, "\n")
  cat("Quantile model:", x$qmodel, "\n")
  cat("Expected shortfall model:", x$emodel, "\n")
  beta <- x$beta.hat
  theta <- x$theta.hat
  names(beta) <- paste0("beta", seq_along(beta) - 1L)
  names(theta) <- paste0("theta", seq_along(theta) - 1L)
  cat("\nQuantile coefficients:\n")
  print(round(beta, digits = digits))
  cat("\nWeighted expected shortfall coefficients:\n")
  print(round(theta, digits = digits))
  invisible(x)
}
