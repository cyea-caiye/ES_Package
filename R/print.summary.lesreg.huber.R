#' Print Adaptive Huber Expected Shortfall Regression Summary
#'
#' @description
#' Print method for objects of class \code{"summary.lesreg.huber"}.
#'
#' @param x an object of class \code{"summary.lesreg.huber"}, usually
#'   produced by \code{\link{summary.lesreg.huber}}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @seealso
#' \code{\link{summary.lesreg.huber}}, \code{\link{lesreg.huber}}
#'
#' @export

print.summary.lesreg.huber <- function(x, digits = max(5L, getOption("digits") - 2L), ...){
  cat("\nCall: ")
  dput(x$call)
  quantile.level <- if(identical(x$tail, "upper")){
    1 - x$tau
  } else{
    x$tau
  }
  cat("\nQuantile level:\n")
  print(format(quantile.level, digits = digits), quote = FALSE)
  cat("\nAlpha3:\n")
  print(format(x$alpha3, digits = digits), quote = FALSE)
  cat("\ngamma:\n")
  print(format(x$gamma, digits = digits), quote = FALSE)
  cat("\nConfidence level:\n")
  print(format(x$level, digits = digits), quote = FALSE)
  cat("\nExpected shortfall regression coefficients:\n")
  coef_print <- x$coefficients
  p_col <- match("P-value", colnames(coef_print))
  if(!is.na(p_col)){
    numeric_part <- coef_print[, -p_col, drop = FALSE]
    numeric_part <- format(round(numeric_part, digits = digits), scientific = FALSE, trim = TRUE)
    p_part <- format.pval(coef_print[, p_col], digits = digits, eps = 10^(-digits))
    coef_print <- cbind(numeric_part, `P-value` = p_part)
  } else{
    coef_print <- format(round(coef_print, digits = digits), trim = TRUE)
  }
  print(coef_print, quote = FALSE)
  invisible(x)
}
