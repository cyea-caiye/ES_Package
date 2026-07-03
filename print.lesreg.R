#' @rdname lesreg
#' @export
print.lesreg <- function(x, ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl)
  }
  cat("\nQuantile regression coefficients:\n")
  print(x$beta.hat, ...)
  cat("\nExpected shortfall regression coefficients:\n")
  print(x$theta.hat, ...)
  nobs <- length(stats::residuals(x))
  rdf <- x$df.residual
  cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
  if(!is.null(na.action <- x$na.action)){
    cat("  (", stats::naprint(na.action), ")\n", sep = "")
  }
  invisible(x)
}
