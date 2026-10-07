#' Summarizing Robust Expected Shortfall Regression Fits
#'
#' @description
#' Summarizes a fitted robust linear expected shortfall (ES) regression model
#' returned by \code{\link{lesreg.huber}}. The method reports coefficient
#' estimates, robust standard errors, confidence intervals, and \eqn{p}-values.
#' The covariance estimator follows the robust inference procedure of
#' He et al. (2023).
#'
#' @param object an object of class \code{"lesreg.huber"}, usually a result
#'   of a call to \code{\link{lesreg.huber}}.
#' @param rob positive numeric scalar used in the default calibration of the
#'   second robustification parameter \code{gamma}. The default is 1.
#' @param gamma positive numeric scalar giving the truncation parameter used
#'   in the robust covariance estimator. If \code{NULL}, it is set to
#'   \code{(rob * n / (p * log(n)))^(1/3)}.
#' @param covariance logical flag indicating whether the estimated covariance
#'   matrix and its sandwich components should be returned.
#' @param level confidence level for the reported confidence intervals. This
#'   must be a single number strictly between 0 and 1.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' A fitted \code{"lesreg.huber"} object contains the first-step quantile
#' regression coefficients and the second-step adaptive Huber ES regression
#' coefficients. This summary method performs statistical inference for the
#' ES regression coefficients.
#'
#' The covariance estimator has a sandwich form. The middle matrix
#' \eqn{\boldsymbol{\Omega}} is estimated robustly by truncating the
#' estimating-equation residuals at the level \code{gamma}, following
#' He et al. (2023). If \code{gamma = NULL}, the truncation level is calibrated
#' automatically as \code{(rob * n / (p * log(n)))^(1/3)}, where \code{n} is
#' the sample size and \code{p} is the number of fitted ES regression
#' coefficients.
#'
#' The reported confidence intervals and \eqn{p}-values are based on a
#' standard normal approximation.
#'
#' @return
#' \code{summary.lesreg.huber} returns an object of class
#' \code{"summary.lesreg.huber"} containing:
#' \describe{
#'   \item{\code{call}}{the call that produced the fitted model.}
#'   \item{\code{terms}}{the \code{\link{terms}} object used in fitting.}
#'   \item{\code{coefficients.results}}{a matrix containing the estimated robust
#'   expected shortfall regression coefficients, their standard errors, lower
#'   and upper confidence limits, and \eqn{p}-values.}
#'   \item{\code{residuals}}{the residuals from the second-step adaptive Huber
#'   regression.}
#'   \item{\code{tau}}{the tail probability used for estimation.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#'   \item{\code{level}}{the confidence level used for the reported confidence
#'   intervals.}
#'   \item{\code{rob}}{the calibration constant used for the default choice of
#'   \code{gamma}.}
#'   \item{\code{gamma}}{the truncation parameter used in the robust covariance
#'   estimator.}
#'   \item{\code{cov}}{the estimated covariance matrix, if
#'   \code{covariance = TRUE}.}
#'   \item{\code{Sigma_hat}}{the estimated outer matrix in the sandwich
#'   covariance calculation, if \code{covariance = TRUE}.}
#'   \item{\code{Omega_gamma_hat}}{the robustly estimated middle matrix in the
#'   sandwich covariance calculation, if \code{covariance = TRUE}.}
#' }
#'
#' @references
#' He, X., Tan, K. M., and Zhou, W.-X. (2023).
#' Robust estimation and inference for expected shortfall regression with many
#' regressors.
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology},
#' 85(4), 1223--1246.
#'
#' @seealso
#' \code{\link{lesreg.huber}},
#' \code{\link{lesreg.huber.fit}},
#' \code{\link{summary.lesreg}}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' tau <- 0.1
#' df <- 4
#'
#' ## Generate data with heavy-tailed errors
#' eps <- rt(n, df = df)
#' y <- 1 + 2 * x1 - x2 + eps
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Fit lower-tail robust ES regression
#' fit <- lesreg.huber(
#'   y ~ x1 + x2,
#'   data = dat,
#'   tau = tau
#' )
#' summary(fit)
#'
#' @export

summary.lesreg.huber <- function(object, rob = 1, gamma = NULL, covariance = FALSE, level = 0.95, ...){
  if(!inherits(object, "lesreg.huber")){
    stop("'object' must be a lesreg.huber object")
  }
  if(!is.numeric(rob) || length(rob) != 1L || is.na(rob) || !is.finite(rob) || rob <= 0){
    stop("'rob' must be a positive numeric scalar")
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
  theta_hat <- object$coefficients
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
    gamma <- (rob * n / (p * log(n)))^(1 / 3)
  }
  quant_residual <- as.vector(Y - X %*% beta_hat)
  indicator <- object$indicator
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
    coefficients.results = coef,
    residuals = object$residuals,
    tau = tau,
    tail = tail,
    level = level,
    rob = rob,
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
