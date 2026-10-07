#' Printing Autoregressive Quantile-Expected Shortfall Fits
#'
#' @description
#' \code{print.estools.ts} prints the model specification and coefficient
#' estimates from a fitted autoregressive quantile-expected shortfall model.
#'
#' @param x an object of class \code{"estools.ts"}.
#' @param digits number of significant digits used for printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The printed output includes the matched function call, quantile level,
#' conditional quantile and expected shortfall model specifications,
#' quantile coefficient estimates, and the final weighted expected shortfall
#' coefficient estimates.
#'
#' @return
#' Invisibly returns \code{x}.
#'
#' @seealso
#' \code{\link{estools.ts}},
#' \code{\link{summary.estools.ts}}
#'
#' @export
print.estools.ts <- function(
    x,
    digits = max(5L, getOption("digits") - 2L),
    ...){
  if(!inherits(x, "estools.ts")){
    stop("'x' must be a estools.ts object")
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\nTau:", x$tau, "\n")
  cat("Quantile model:", x$qmodel, "\n")
  cat("Expected shortfall model:", x$emodel, "\n")
  beta <- x$beta.hat
  theta <- x$coefficients
  names(beta) <- paste0("beta", seq_along(beta) - 1L)
  names(theta) <- paste0("theta", seq_along(theta) - 1L)
  cat("\nQuantile coefficients:\n")
  print(round(beta, digits = digits))
  cat("\nWeighted expected shortfall coefficients:\n")
  print(round(theta, digits = digits))
  invisible(x)
}
