#' @rdname lesreg.highdim
#' @export
print.summary.lesreg.highdim <- function(x, ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
  }
  cat("\nHigh-dimensional expected shortfall regression\n")
  if(!is.null(x$tau)){
    cat("\nTau:", x$tau, "\n")
  }
  if(!is.null(x$uppertail)){
    cat("Tail:", ifelse(x$uppertail, "upper", "lower"), "\n")
  }
  if(!is.null(x$tuning.method)){
    cat("Tuning method:", x$tuning.method, "\n")
  }
  if(!is.null(x$standardize)){
    cat("Standardize:", x$standardize, "\n")
  }
  cat("\nTuning parameters:\n")
  cat("  Quantile regression lambda:", x$lambda_beta, "\n")
  cat("  Expected shortfall regression lambda:", x$lambda_theta, "\n")
  cat("\nNumber of nonzero coefficients:\n")
  cat("  Quantile regression:", x$nzero_beta, "of", x$p, "\n")
  cat("  Expected shortfall regression:", x$nzero_theta, "of", x$p, "\n")
  cat("\nQuantile regression coefficients:\n")
  print(x$quantile.coefficients[x$quantile.coefficients != 0], ...)
  cat("\nExpected shortfall regression coefficients:\n")
  print(x$es.coefficients[x$es.coefficients != 0], ...)
  invisible(x)
}
