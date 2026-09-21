#' Summarizing Autoregressive Quantile-Expected Shortfall Fits
#'
#' @description
#' Summary method for objects of class \code{"lesreg.ts"}.
#'
#' @param object an object of class \code{"lesreg.ts"}.
#' @param level confidence level for the reported confidence intervals.
#' @param covariance logical. If \code{TRUE}, the estimated covariance
#'   matrix is included in the returned object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.lesreg.ts"}.
#'
#' @export
summary.lesreg.ts <- function(object,
                              level = 0.95,
                              covariance = FALSE,
                              ...){

  if(!inherits(object, "lesreg.ts")){
    stop("'object' must be a lesreg.ts object")
  }
  if(!is.numeric(level) || length(level) != 1L || is.na(level) || level <= 0 || level >= 1){
    stop("'level' must be a single number between 0 and 1")
  }

  if(!is.logical(covariance) || length(covariance) != 1L || is.na(covariance)){
    stop("'covariance' must be TRUE or FALSE")
  }
  theta_hat <- object$theta.hat
  VC <- object$VC
  n <- length(object$quantile.fitted)
  robust_se <- sqrt(diag(VC) / n)
  z_alpha <- stats::qnorm(1 - (1 - level) / 2)
  lower <- theta_hat - z_alpha * robust_se
  upper <- theta_hat + z_alpha * robust_se
  pivot <- theta_hat / robust_se
  p_values <- 2 * stats::pnorm(abs(pivot), lower.tail = FALSE)
  coef_names <- names(theta_hat)
  if(is.null(coef_names)){
    coef_names <- paste0("theta", seq_along(theta_hat) - 1L)
  }
  coefficients <- cbind(
    Estimate = theta_hat,
    `Std. Error` = robust_se,
    `Lower CI` = lower,
    `Upper CI` = upper,
    `P-value` = p_values
  )
  rownames(coefficients) <- coef_names
  ans <- list(
    call = object$call,
    coefficients = coefficients,
    beta.hat = object$beta.hat,
    theta.unweighted = object$theta.unweighted,
    tau = object$tau,
    qmodel = object$qmodel,
    emodel = object$emodel,
    level = level,
    nobs = n
  )

  if(covariance){
    ans$cov <- VC / n
  }

  class(ans) <- "summary.lesreg.ts"
  ans
}
