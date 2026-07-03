#' Summarizing Linear Expected Shortfall Regression Fits
#'
#' @description
#' Summary method for class \code{"lesreg"}.
#'
#' @usage
#' \method{summary}{lesreg}(object,
#'         method = c("default", "ind", "scl_N", "scl_sp", "rwb"),
#'         covariance = FALSE, level = 0.95, B = 500, ...)
#'
#' \method{print}{summary.lesreg}(x,
#'         digits = max(5L, getOption("digits") - 2L), ...)
#'
#' @param object an object of class \code{"lesreg"}, usually a result of a call
#'   to \code{\link{lesreg}}.
#' @param x an object of class \code{"summary.lesreg"}, usually a result of a
#'   call to \code{summary.lesreg}.
#' @param method character string specifying the method used to compute
#'   standard errors. There are currently five available methods:
#'   \describe{
#'     \item{\code{"default"}}{the default sandwich estimator for the two-step
#'     expected shortfall regression fit.}
#'     \item{\code{"ind"}}{an estimator based on a conditional truncated
#'     variance approximation which is treated as independent of the
#'     covariates.}
#'     \item{\code{"scl_N"}}{an estimator based on a normal scale
#'     approximation to the conditional truncated variance.}
#'     \item{\code{"scl_sp"}}{an estimator based on a semiparametric scale
#'     approximation to the conditional truncated variance.}
#'     \item{\code{"rwb"}}{a random-weighting bootstrap estimator using
#'     exponential random weights.}
#'   }
#' @param covariance logical flag indicating whether the estimated covariance
#'   matrix and, where available, its sandwich components should be returned.
#' @param level confidence level for the reported confidence intervals.
#' @param B number of random-weighting bootstrap replications. Used only when
#'   \code{method = "rwb"}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The fitted \code{"lesreg"} object contains both the first-step quantile
#' regression coefficients and the second-step expected shortfall regression
#' coefficients. The summary method reports inference for the expected
#' shortfall regression coefficients.
#'
#' When \code{method = "default"}, standard errors are computed from the
#' sandwich form used for the two-step expected shortfall regression estimator.
#' The methods \code{"ind"}, \code{"scl_N"} and \code{"scl_sp"} use alternative
#' estimates of the conditional truncated variance term. These estimates are
#' computed by \code{\link[esreg]{conditional_truncated_variance}}.
#'
#' When \code{method = "rwb"}, standard errors and confidence intervals are
#' computed from random-weighting bootstrap replications. The argument
#' \code{B} controls the number of bootstrap replications.
#'
#' P-values are computed using the standard normal approximation.
#'
#' @return
#' The function \code{summary.lesreg} computes and returns a list of summary
#' statistics of the fitted expected shortfall regression model given in
#' \code{object}, using the components \code{call}, \code{terms},
#' \code{x.fit}, \code{y.fit}, \code{beta.hat} and \code{theta.hat} from its
#' argument, plus
#' \describe{
#'   \item{coefficients}{a matrix with columns for the estimated coefficients,
#'   their standard errors, lower and upper confidence limits, and p-values.}
#'   \item{residuals}{the residuals from the fitted second-step regression.}
#'   \item{rdf}{the residual degrees of freedom.}
#'   \item{tau}{the tail probability used for estimation.}
#'   \item{method}{the variance estimation method used.}
#'   \item{tail}{a character string indicating whether the lower or upper tail was fitted.}
#'   \item{level}{the confidence level used for the reported confidence intervals.}
#'   \item{cov}{the estimated covariance matrix, if \code{covariance = TRUE}.}
#'   \item{Sigma_hat}{the estimated outer matrix in the sandwich covariance
#'   calculation, if \code{covariance = TRUE} and \code{method != "rwb"}.}
#'   \item{Omega_hat}{the estimated middle matrix in the sandwich covariance
#'   calculation, if \code{covariance = TRUE} and \code{method != "rwb"}.}
#' }
#'
#' @seealso
#' \code{\link{lesreg}}, \code{\link{lesreg.fit}}, \code{\link{lesreg.wfit}},
#' \code{\link[esreg]{conditional_truncated_variance}}
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05)
#' summary(fit)
#' summary(fit, method = "ind")
#' summary(fit, method = "scl_N", covariance = TRUE)
#' summary(fit, method = "rwb", B = 100)
#'
#' @export


summary.lesreg <- function(object, method = c("default", "ind", "scl_N", "scl_sp", "rwb"), covariance = FALSE, level = 0.95, B = 500, ...){
  if(!inherits(object, "lesreg"))
    stop("'object' must be a lesreg object")
  if(!is.logical(covariance) || length(covariance) != 1L || is.na(covariance))
    stop("'covariance' must be TRUE or FALSE")
  if(!is.numeric(level) || length(level) != 1L || is.na(level) || level <= 0 || level >= 1)
    stop("'level' must be a single number between 0 and 1")
  method <- match.arg(method)
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
  rdf <- object$df.residual
  if(p == 0L){
    stop("empty model has no coefficients to summarize")
  }
  vnames <- names(theta_hat)
  if(is.null(vnames)) vnames <- colnames(X)
  if(is.null(vnames)) vnames <- paste0("x", seq_len(p))
  I <- object$indicator
  quant_residual <- as.vector(Y - X %*% beta_hat)
  Sigma_hat <- (1 / n) * t(X) %*% X
  if(method == "default"){
    omega_hat <- quant_residual * I + tau * as.vector(X %*% (beta_hat - theta_hat))
    Omega_hat <- (1 / n) * t(X) %*% (X * as.vector(omega_hat^2))
    sandwich_var <- solve(Sigma_hat) %*% Omega_hat %*% solve(Sigma_hat)
    robust_se <- sqrt(diag(sandwich_var)) / (tau * sqrt(n))
    pivot <- theta_hat / robust_se
    z_alpha <- stats::qnorm(1 - (1 - level) / 2)
    lower <- theta_hat - z_alpha * robust_se
    upper <- theta_hat + z_alpha * robust_se
  } else if(method == "rwb"){
    if(!is.numeric(B) || length(B) != 1L || is.na(B) || B <= 1 || B != as.integer(B)){
      stop("'B' must be an integer greater than 1")
    }
    coef.rwb <- matrix(NA, nrow = B, ncol = p)
    data.weights <- object$weights
    if(is.null(data.weights)){
      data.weights <- rep(1, n)
    }
    for(b in 1:B){
      e <- stats::rexp(n, rate = 1)
      rweights <- data.weights * e
      fit.b <- lesreg.wfit(x = X, y = Y,weights = rweights, tau = tau, method = object$method, upperTail = FALSE,
                           eps = object$eps, lm.tol = object$lm.tol, singular.ok = object$singular.ok, ...)
      coef.rwb[b, ] <- fit.b$theta.hat
    }
    diff <- sqrt(n) * sweep(coef.rwb, 2, theta_hat, FUN = "-")
    sandwich_var <- stats::cov(diff) / n
    robust_se <- sqrt(diag(sandwich_var))
    pivot <- theta_hat / robust_se
    q.low <- apply(diff, 2, stats::quantile, probs = (1 - level) / 2, names = FALSE)
    q.high <- apply(diff, 2, stats::quantile, probs = 1 - (1 - level) / 2, names = FALSE)
    lower <- theta_hat - q.high / sqrt(n)
    upper <- theta_hat - q.low / sqrt(n)
  } else{
    V <- esreg::conditional_truncated_variance(y = quant_residual, x = X, approach = method)
    conditional_hetero <- V / tau + (1 - tau) / tau * as.vector(X %*% beta_hat - X %*% theta_hat)^2
    Omega_hat <- (1 / n) * t(X) %*% (X * as.vector(conditional_hetero))
    sandwich_var <- solve(Sigma_hat) %*% Omega_hat %*% solve(Sigma_hat)
    robust_se <- sqrt(diag(sandwich_var)) / sqrt(n)
    pivot <- theta_hat / robust_se
    z_alpha <- stats::qnorm(1 - (1 - level) / 2)
    lower <- theta_hat - z_alpha * robust_se
    upper <- theta_hat + z_alpha * robust_se
  }
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
    rdf = rdf,
    tau = tau,
    method = method,
    tail = tail,
    level = level
  )
  if(covariance){
    ans$cov <- sandwich_var
    if(method != "rwb"){
      ans$Sigma_hat <- Sigma_hat
      ans$Omega_hat <- Omega_hat
    }
  }
  class(ans) <- "summary.lesreg"
  ans
}
