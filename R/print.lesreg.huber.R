#' @rdname lesreg.huber
#' @export
print.lesreg.huber <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  cat(
    "\nCall:\n",
    paste(
      deparse(x$call),
      sep = "\n",
      collapse = "\n"
    ),
    "\n\n",
    sep = ""
  )
  quantile.level <- if(identical(x$tail, "upper")){
    1 - x$tau
  } else{
    x$tau
  }
  cat("Quantile level:\n")
  print.default(
    format(quantile.level, digits = digits),
    quote = FALSE
  )
  cat("\n")
  if(length(x$beta.hat)){
    cat("Quantile regression coefficients:\n")
    print.default(
      format(x$beta.hat, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  } else{
    cat("No quantile regression coefficients\n")
  }
  cat("\n")
  if(length(x$theta.hat)){
    cat("Expected shortfall regression coefficients:\n")
    print.default(
      format(x$theta.hat, digits = digits),
      print.gap = 2L,
      quote = FALSE
    )
  } else{
    cat("No expected shortfall regression coefficients\n")
  }
  cat("\n")
  cat("Adaptive Huber threshold:\n")
  print.default(
    format(x$huber.tau, digits = digits),
    quote = FALSE
  )
  cat("\n")
  cat("Adaptive Huber iterations:\n")
  print.default(
    x$huber.iteration,
    quote = FALSE
  )
  cat("\n")
  invisible(x)
}
