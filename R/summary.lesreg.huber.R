#' Summarizing Adaptive Huber Expected Shortfall Regression Fits
#'
#' @description
#' Summary method for objects of class \code{"lesreg.huber"}.
#'
#' @param object an object of class \code{"lesreg.huber"}, usually a result
#'   of a call to \code{\link{lesreg.huber}}.
#' @param alpha3 positive numeric scalar used in the default calibration of
#'   the second robustification parameter \code{gamma}. The default is 1.
#' @param gamma positive numeric scalar giving the truncation parameter used
#'   in robust covariance estimation. If \code{NULL}, it is set to
#'   \code{(alpha3 * n / (p * log(n)))^(1/3)}.
#' @param covariance logical flag indicating whether the estimated covariance
#'   matrix and its sandwich components should be returned.
#' @param level confidence level for the reported confidence intervals.
#' @param ... further arguments passed to or from other methods.
#'
#' @return
#' A list of class \code{"summary.lesreg.huber"} containing coefficient
#' estimates, standard errors, confidence intervals, p-values, and the
#' robust covariance tuning parameters.
#'
#' @export

summary.lesreg.huber <- function(object, alpha3 = 1, gamma = NULL, covariance = FALSE, level = 0.95, ...){
  if(!inherits(object, "lesreg.huber")){
    stop("'object' must be a lesreg.huber object")
  }
  if(!is.numeric(alpha3) || length(alpha3) != 1L || is.na(alpha3) || !is.finite(alpha3) || alpha3 <= 0){
    stop("'alpha3' must be a positive numeric scalar")
  }
  if(!is.null(gamma)){
    if(!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || !is.finite(gamma) || gamma <= 0){
      stop("'gamma' must be NULL or a positive numeric scalar")
    }
  }
  if(!is.logical(covariance) || length(covariance) != 1L || is.na(covariance)){
    stop("'covariance' must be TRUE or FALSE")
  }
  if(!is.numeric(level) || length(level) != 1L || is.na(level) || !is.finite(level) || level <= 0 || level >= 1){
    stop("'level' must be a single number between 0 and 1")
  }
  call <- object$call
  terms <- object$terms
  tail <- object$tail
  X <- object$x.fit
  Y <- object$y.fit
  tau <- object$tau
  beta_hat <- object$beta.hat
  theta_hat <- object$theta.hat
  n <- nrow(X)
  p <- length(theta_hat)
  if(p == 0L){
    stop("empty model has no coefficients to summarize")
  }
  vnames <- names(theta_hat)
  if(is.null(vnames)){
    vnames <- colnames(X)
  }
  if(is.null(vnames)){
    vnames <- paste0("x", seq_len(p))
  }
  if(is.null(gamma)){
    gamma <- (alpha3 * n / (p * log(n)))^(1 / 3)
  }
  quant_residual <- as.vector(Y - X %*% beta_hat)
  indicator <- as.numeric(quant_residual <= 0)
  omega_hat <- quant_residual * indicator + tau * as.vector(X %*% (beta_hat - theta_hat))
  psi_gamma <- sign(omega_hat) * pmin(abs(omega_hat), gamma)
  Sigma_hat <- crossprod(X) / n
  Omega_gamma_hat <- crossprod(X, X * as.vector(psi_gamma^2)) / n
  Sigma_inv <- solve(Sigma_hat)
  sandwich_core <- Sigma_inv %*% Omega_gamma_hat %*% Sigma_inv
  cov_hat <- sandwich_core / (tau^2 * n)
  robust_se <- sqrt(diag(cov_hat))
  pivot <- theta_hat / robust_se
  z_alpha <- stats::qnorm(1 - (1 - level) / 2)
  lower <- theta_hat - z_alpha * robust_se
  upper <- theta_hat + z_alpha * robust_se
  p_values <- 2 * stats::pnorm(abs(pivot), lower.tail = FALSE)
  coef <- cbind(
    Estimate = theta_hat,
    `Std. Error` = robust_se,
    `Lower CI` = lower,
    `Upper CI` = upper,
    `P-value` = p_values
  )
  dimnames(coef) <- list(vnames, colnames(coef))
  ans <- list(
    call = call,
    terms = terms,
    coefficients = coef,
    residuals = object$residuals,
    tau = tau,
    tail = tail,
    level = level,
    alpha3 = alpha3,
    gamma = gamma
  )
  if(covariance){
    ans$cov <- cov_hat
    ans$Sigma_hat <- Sigma_hat
    ans$Omega_gamma_hat <- Omega_gamma_hat
  }
  class(ans) <- "summary.lesreg.huber"
  ans
}
