#' Summarizing Autoregressive Quantile-Expected Shortfall Fits
#'
#' @description
#' Summary method for objects of class \code{"lesreg.ts"}.
#'
#' @details
#' Inference is reported for the final weighted expected shortfall estimator.
#' Standard errors are based on the sandwich covariance matrix computed during
#' model fitting using the selected \code{covMethod} and
#' \code{truncationMethod}.
#'
#' @param object an object of class \code{"lesreg.ts"}.
#' @param level confidence level for the reported confidence intervals.
#' @param covariance logical. If \code{TRUE}, the estimated covariance matrix
#'   together with the sandwich components \code{Sigma_hat} and
#'   \code{Omega_hat} is included in the returned object.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \code{"summary.lesreg.ts"} containing the
#'   weighted expected shortfall coefficient estimates, standard errors,
#'   confidence intervals, and p-values. If \code{covariance = TRUE}, the
#'   estimated covariance matrix and sandwich components are also returned.
#'
#' @examples
#' set.seed(123)
#' y <- as.numeric(arima.sim(model = list(ar = 0.3), n = 500))
#'
#' fit <- lesreg.ts(
#'   y,
#'   tau = 0.1,
#'   qmodel = "SAV",
#'   emodel = "SAV"
#' )
#'
#' summary(fit)
#' summary(fit, covariance = TRUE)
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
    tau = object$tau,
    qmodel = object$qmodel,
    emodel = object$emodel,
    covMethod = object$covMethod,
    truncationMethod = object$truncationMethod,
    level = level,
    nobs = n
  )

  if(covariance){
    ans$cov <- VC / n
    ans$Sigma_hat <- object$weighted.fit$Sigma
    ans$Omega_hat <- object$weighted.fit$Omega
  }

  class(ans) <- "summary.lesreg.ts"
  ans
}
