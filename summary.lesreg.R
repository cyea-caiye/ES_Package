#' Summary Method for Linear Expected Shortfall Regression
#'
#' @description
#' Returns a summary list for a fitted linear expected shortfall regression
#' model.
#'
#' The function computes standard errors, confidence intervals, and p-values
#' for the expected shortfall regression coefficients using one of several
#' sandwich-type variance estimation methods.
#'
#' @usage
#' \method{summary}{lesreg}(object, method = c("default", "ind", "scl_N", "scl_sp"), covariance = FALSE, level = 0.95, ...)
#'
#' @param object An object of class \code{"lesreg"}, typically produced by a
#'   call to \code{\link{lesreg}}.
#' @param method Character string specifying the method used to compute
#'   standard errors. The available methods are:
#'   \describe{
#'     \item{\code{"default"}}{The default sandwich estimator based on the
#'     two-step expected shortfall regression framework. It uses the empirical
#'     analogue of \eqn{\Omega = E(\omega^2XX^\top)}, where
#'     \eqn{\omega = (Y-X^\top\beta)1(Y \leq X^\top\beta) + \tau X^\top(\beta-\theta)}.}
#'     \item{\code{"ind"}}{An estimator based on the joint quantile and
#'     expected shortfall regression framework. The conditional truncated
#'     variance term is estimated by the sample variance of the negative
#'     quantile residuals and is treated as independent of the covariates.}
#'     \item{\code{"scl_N"}}{An estimator based on the joint quantile and
#'     expected shortfall regression framework. The conditional truncated
#'     variance is estimated by a normal scale approximation after estimating
#'     the conditional location and scale of the quantile residuals.}
#'     \item{\code{"scl_sp"}}{An estimator based on the joint quantile and
#'     expected shortfall regression framework. The conditional truncated
#'     variance is estimated by a semiparametric scale approach using the
#'     standardized quantile residual distribution.}
#'   }
#' @param covariance Logical flag. If \code{TRUE}, the matrix
#'   \code{sandwich_var} and its components are returned. The default is
#'   \code{FALSE}.
#' @param level Confidence level for the reported Wald-type confidence
#'   intervals. The default is \code{0.95}.
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' The fitted \code{"lesreg"} object contains both the first-step quantile
#' regression coefficients and the second-step expected shortfall regression
#' coefficients. This summary method focuses on inference for the expected
#' shortfall coefficients.
#'
#' The \code{"default"} method follows the two-step expected shortfall
#' regression inference based on the sandwich form involving
#' \eqn{\Sigma = E(XX^\top)} and
#' \eqn{\Omega = E(\omega^2XX^\top)}.
#'
#' The methods \code{"ind"}, \code{"scl_N"}, and \code{"scl_sp"} follow the
#' joint quantile and expected shortfall regression covariance structure. These
#' methods differ in how they estimate the conditional truncated variance term
#' \deqn{
#' Var(Y-X^\top\beta \mid Y \leq X^\top\beta, X).
#' }
#' The actual computation of this quantity is delegated to
#' \code{\link[esreg]{conditional_truncated_variance}}.
#'
#' P-values are computed from the standard normal approximation.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{terms}{The \code{\link{terms}} object used in the fitted model.}
#'   \item{coefficients}{A matrix containing the estimated expected shortfall
#'   regression coefficients, their estimated standard errors, confidence
#'   intervals, and associated p-values.}
#'   \item{residuals}{Residuals from the fitted expected shortfall regression
#'   model.}
#'   \item{rdf}{Residual degrees of freedom.}
#'   \item{tau}{The tail probability used for estimation.}
#'   \item{method}{The variance estimation method used.}
#'   \item{tail}{Character string indicating whether the fitted model
#'   corresponds to the lower or upper tail.}
#'   \item{level}{The confidence level used for the reported confidence
#'   intervals.}
#'   \item{cov}{The matrix \code{sandwich_var} returned by the selected
#'   variance estimation method. This component is returned only when
#'   \code{covariance = TRUE}.}
#'   \item{Sigma_hat}{The estimated matrix \eqn{n^{-1}X^\top X}. This
#'   component is returned only when \code{covariance = TRUE}.}
#'   \item{Omega_hat}{The estimated middle matrix in the sandwich variance
#'   calculation. This component is returned only when
#'   \code{covariance = TRUE}.}
#' }
#'
#' @seealso
#' \code{\link{lesreg}}, \code{\link{lesreg.fit}},
#' \code{\link[quantreg]{rq.fit}},
#' \code{\link[esreg]{conditional_truncated_variance}}
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' fit <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05)
#'
#' summary(fit)
#' summary(fit, method = "ind")
#' summary(fit, method = "scl_N", covariance = TRUE)
#'
#' @export

summary.lesreg <- function(object, method = c("default", "ind", "scl_N", "scl_sp"), covariance = FALSE, level = 0.95, ...){
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
  } else{
    V <- esreg::conditional_truncated_variance(y = quant_residual, x = X, approach = method)
    conditional_hetero <- V / tau + (1 - tau) / tau * as.vector(X %*% beta_hat - X %*% theta_hat)^2
    Omega_hat <- (1 / n) * t(X) %*% (X * as.vector(conditional_hetero))
    sandwich_var <- solve(Sigma_hat) %*% Omega_hat %*% solve(Sigma_hat)
    robust_se <- sqrt(diag(sandwich_var)) / sqrt(n)
    pivot <- theta_hat / robust_se
  }
  z_alpha <- stats::qnorm(1 - (1 - level) / 2)
  p_values <- 2 * stats::pnorm(abs(pivot), lower.tail = FALSE)
  lower <- theta_hat - z_alpha * robust_se
  upper <- theta_hat + z_alpha * robust_se
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
    ans$Sigma_hat <- Sigma_hat
    ans$Omega_hat <- Omega_hat
  }
  class(ans) <- "summary.lesreg"
  ans
}
