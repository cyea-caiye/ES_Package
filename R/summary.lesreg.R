#' Summarizing Linear Expected Shortfall Regression Fits
#'
#' @description
#' Summarizes a fitted linear expected shortfall (ES) regression model and
#' provides standard errors, confidence intervals, and \eqn{p}-values for the
#' estimated ES regression coefficients.
#'
#' @param x an object of class \code{"summary.lesreg"}, usually a result of a
#'   call to \code{summary.lesreg}.
#' @param object an object of class \code{"lesreg"}, usually a result of a call
#'   to \code{\link{lesreg}}.
#' @param method character string specifying the method used to estimate the
#'   covariance matrix of the ES regression coefficients. There are currently
#'   five available methods:
#'   \describe{
#'     \item{\code{"default"}}{the default sandwich covariance estimator for
#'     the two-step expected shortfall regression estimator.}
#'     \item{\code{"ind"}}{an estimator that treats the conditional truncated
#'     variance as constant across covariate values.}
#'     \item{\code{"scl_N"}}{an estimator based on a location-scale model with
#'     a normal approximation to the conditional truncated variance.}
#'     \item{\code{"scl_sp"}}{an estimator based on a location-scale model with
#'     a semiparametric approximation to the conditional truncated variance.}
#'     \item{\code{"rwb"}}{a random-weighting bootstrap estimator using
#'     independent exponential random weights.}
#'   }
#' @param covariance logical flag indicating whether the estimated covariance
#'   matrix and, where available, its sandwich components should be returned.
#' @param level confidence level for the reported confidence intervals. This
#'   must be a single number strictly between 0 and 1.
#' @param B number of random-weighting bootstrap replications. This must be an
#'   integer greater than 1 and is used only when \code{method = "rwb"}.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' A fitted \code{"lesreg"} object contains both the first-step quantile
#' regression coefficients and the second-step expected shortfall regression
#' coefficients. The summary method performs statistical inference for the
#' expected shortfall regression coefficients. See \code{\link{lesreg}} for
#' details on the underlying two-step estimation procedure.
#'
#' When \code{method = "default"}, standard errors are computed using the
#' sandwich covariance estimator associated with the two-step expected
#' shortfall regression estimator.
#'
#' The methods \code{"ind"}, \code{"scl_N"}, and \code{"scl_sp"} use
#' alternative estimators of the conditional truncated variance. The
#' \code{"ind"} method treats this variance as constant across covariate
#' values, whereas \code{"scl_N"} and \code{"scl_sp"} use location-scale
#' approximations. These approaches follow the conditional truncated variance
#' estimation methods described by Dimitriadis and Bayer (2019).
#'
#' When \code{method = "rwb"}, the covariance matrix is estimated using a
#' random-weighting bootstrap with independent exponential random weights.
#' The reported standard errors are obtained from the bootstrap covariance
#' matrix, and confidence intervals are constructed using the empirical
#' quantiles of the centered bootstrap coefficient estimates. The argument
#' \code{B} controls the number of bootstrap replications.
#'
#' For all methods, the reported \eqn{p}-values are computed using a standard
#' normal approximation.
#'
#' @return
#' \code{summary.lesreg} returns an object of class
#' \code{"summary.lesreg"} containing summary statistics and inference results
#' for the fitted expected shortfall regression model. The documented
#' components are:
#' \describe{
#'   \item{\code{coefficients.results}}{a matrix containing the estimated
#'   expected shortfall regression coefficients, their standard errors, lower
#'   and upper confidence limits, and \eqn{p}-values.}
#'   \item{\code{residuals}}{the residuals from the second-step least-squares
#'   regression, defined as the pseudo-response minus the fitted expected
#'   shortfall values.}
#'   \item{\code{rdf}}{the residual degrees of freedom.}
#'   \item{\code{tau}}{the tail probability used for estimation.}
#'   \item{\code{method}}{the method used to estimate the covariance matrix.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#'   \item{\code{level}}{the confidence level used for the reported confidence
#'   intervals.}
#'   \item{\code{cov}}{the estimated covariance matrix, if
#'   \code{covariance = TRUE}.}
#'   \item{\code{Sigma_hat}}{the estimated outer matrix in the sandwich
#'   covariance calculation, if \code{covariance = TRUE} and
#'   \code{method != "rwb"}.}
#'   \item{\code{Omega_hat}}{the estimated middle matrix in the sandwich
#'   covariance calculation, if \code{covariance = TRUE} and
#'   \code{method != "rwb"}.}
#' }
#'
#' @references
#' Dimitriadis, T., and Bayer, S. (2019).
#' A joint quantile and expected shortfall regression framework.
#' \emph{Electronic Journal of Statistics}, 13(1), 1823--1871.
#'
#' @seealso
#' \code{\link{lesreg}}, \code{\link{lesreg.fit}}, \code{\link{lesreg.wfit}}
#'
#' @examples
#' set.seed(2026)
#' n <- 500
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05)
#'
#' ## The default covariance estimation method.
#' summary(fit)
#'
#' ## The "ind" estimation method in Dimitriadis and Bayer (2019).
#' summary(fit, method = "ind")
#'
#' ## The "scl_N" estimation method in Dimitriadis and Bayer (2019).
#' summary(fit, method = "scl_N", covariance = TRUE)
#'
#' ## The random-weighting bootstrap method.
#' summary(fit, method = "rwb", B = 100)
#'
#' @export


summary.lesreg <- function(
    object,
    method = c("default", "ind", "scl_N", "scl_sp", "rwb"),
    covariance = FALSE,
    level = 0.95,
    B = 500,
    ...
) {
  if (!inherits(object, "lesreg")) {
    stop("'object' must be a lesreg object")
  }

  if (!is.logical(covariance) ||
      length(covariance) != 1L ||
      is.na(covariance)) {
    stop("'covariance' must be TRUE or FALSE")
  }

  if (!is.numeric(level) ||
      length(level) != 1L ||
      is.na(level) ||
      level <= 0 ||
      level >= 1) {
    stop("'level' must be a single number between 0 and 1")
  }

  method <- match.arg(method)

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
  rdf <- object$df.residual

  weights <- if (is.null(object$weights)) {
    rep(1, n)
  } else {
    object$weights
  }

  if (p == 0L) {
    stop("empty model has no coefficients to summarize")
  }

  vnames <- names(theta_hat)
  if (is.null(vnames)) {
    vnames <- colnames(X)
  }
  if (is.null(vnames)) {
    vnames <- paste0("x", seq_len(p))
  }

  if (method == "default") {

    I <- object$indicator
    quant_residual <- as.vector(Y - X %*% beta_hat)

    omega_hat <-
      quant_residual * I +
      tau * as.vector(X %*% (beta_hat - theta_hat))

    Sigma_hat <- crossprod(
      X,
      X * weights
    ) / n

    Omega_hat <- crossprod(
      X,
      X * (weights * omega_hat)^2
    ) / n

    asymptotic_var <-
      solve(Sigma_hat) %*%
      Omega_hat %*%
      solve(Sigma_hat)

    cov_hat <- asymptotic_var / (tau^2 * n)

    robust_se <- sqrt(diag(cov_hat))
    pivot <- theta_hat / robust_se

    z_alpha <- stats::qnorm(1 - (1 - level) / 2)

    lower <- theta_hat - z_alpha * robust_se
    upper <- theta_hat + z_alpha * robust_se

  } else if (method == "rwb") {

    if (!is.numeric(B) ||
        length(B) != 1L ||
        is.na(B) ||
        B <= 1 ||
        B != as.integer(B)) {
      stop("'B' must be an integer greater than 1")
    }

    coef.rwb <- matrix(
      NA_real_,
      nrow = B,
      ncol = p,
      dimnames = list(NULL, vnames)
    )

    for (b in seq_len(B)) {

      e <- stats::rexp(n, rate = 1)
      rweights <- weights * e

      fit.b <- lesreg.wfit(
        x = X,
        y = Y,
        weights = rweights,
        tau = tau,
        method = object$method,
        upperTail = FALSE,
        eps = object$eps,
        lm.tol = object$lm.tol,
        singular.ok = object$singular.ok,
        ...
      )

      coef.rwb[b, ] <- fit.b$coefficients
    }

    diff <- sqrt(n) *
      sweep(coef.rwb, 2, theta_hat, FUN = "-")

    cov_hat <- stats::cov(diff) / n
    robust_se <- sqrt(diag(cov_hat))
    pivot <- theta_hat / robust_se

    q.low <- apply(
      diff,
      2,
      stats::quantile,
      probs = (1 - level) / 2,
      names = FALSE
    )

    q.high <- apply(
      diff,
      2,
      stats::quantile,
      probs = 1 - (1 - level) / 2,
      names = FALSE
    )

    lower <- theta_hat - q.high / sqrt(n)
    upper <- theta_hat - q.low / sqrt(n)

  } else {

    quant_residual <- as.vector(
      Y - X %*% beta_hat
    )

    Sigma_hat <- crossprod(
      X,
      X * weights
    ) / n

    V <- .lesreg_conditional_truncated_variance(
      y = quant_residual,
      x = X,
      approach = method
    )

    conditional_hetero <-
      (
        V / tau +
          (1 - tau) / tau *
          as.vector(
            X %*% beta_hat -
              X %*% theta_hat
          )^2
      ) *
      weights^2

    Omega_hat <- crossprod(
      X,
      X * conditional_hetero
    ) / n

    asymptotic_var <-
      solve(Sigma_hat) %*%
      Omega_hat %*%
      solve(Sigma_hat)

    cov_hat <- asymptotic_var / n

    robust_se <- sqrt(diag(cov_hat))
    pivot <- theta_hat / robust_se

    z_alpha <- stats::qnorm(
      1 - (1 - level) / 2
    )

    lower <- theta_hat - z_alpha * robust_se
    upper <- theta_hat + z_alpha * robust_se
  }

  p_values <- 2 *
    stats::pnorm(
      abs(pivot),
      lower.tail = FALSE
    )

  coef <- cbind(
    Estimate = theta_hat,
    `Std. Error` = robust_se,
    `Lower CI` = lower,
    `Upper CI` = upper,
    `P-value` = p_values
  )

  dimnames(coef) <- list(
    vnames,
    colnames(coef)
  )

  ans <- list(
    call = call,
    terms = terms,
    coefficients.results = coef,
    residuals = object$residuals,
    rdf = rdf,
    tau = tau,
    method = method,
    tail = tail,
    level = level
  )

  if (covariance) {
    ans$cov <- cov_hat

    if (method != "rwb") {
      ans$Sigma_hat <- Sigma_hat
      ans$Omega_hat <- Omega_hat
    }
  }

  class(ans) <- "summary.lesreg"
  ans
}
