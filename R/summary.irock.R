#' Summarize an i-Rock fit
#'
#' @description
#' Summarizes a fitted linear expected shortfall (ES) regression model and
#' provides standard errors, confidence intervals, and \eqn{p}-values for the
#' estimated ES regression coefficients, under the i-Rock framework.
#'
#' @param object A fitted object returned by \code{\link{irock}}.
#' @param x An object returned by \code{summary.irock()}.
#' @param method Covariance estimation method. Only \code{"rwb"} is supported.
#' @param covariance Logical; return the bootstrap covariance matrix?
#' @param level A confidence level strictly between zero and one.
#' @param B Number of bootstrap replications, an integer greater than one.
#' @param digits Number of significant digits when printing the summary.
#' @param ... Reserved in \code{summary.irock()}; arguments passed to
#'   \code{print()} for the formatted coefficient table in the print method.
#'
#' @details
#' Each replication draws independent mean-one exponential weights for the
#' observations retained by \code{irock()}. It refits the initial quantile
#' and expected shortfall estimators and the final integrated quantile loss,
#' holding the fitted covariate bins or nearest-neighbor neighborhoods fixed.
#' The same observation weight is used at every stage and in every
#' neighborhood that contains that observation.
#'
#' Standard errors are bootstrap standard deviations. Confidence intervals
#' are basic bootstrap intervals, formed by subtracting the reversed quantiles
#' of the bootstrap coefficient differences from the original estimates.
#' Two-sided \eqn{p}-values use a normal approximation and are missing for
#' coefficients with zero bootstrap standard error.
#'
#' This implementation uses independent observation weights and is intended
#' for independent observations. It does not adjust for serial or cluster
#' dependence. Set the R random seed before calling \code{summary()} for
#' reproducibility.
#'
#' @return
#' \code{summary.irock} returns an object of class
#' \code{"summary.irock"} containing summary statistics and inference results
#' for the fitted expected shortfall regression model. The documented
#' components are:
#' \describe{
#'   \item{\code{call}}{the original model call.}
#'   \item{\code{terms}}{the terms object for the fitted model.}
#'   \item{\code{coefficients.results}}{a matrix containing the estimated
#'   expected shortfall regression coefficients, their bootstrap standard
#'   errors, lower and upper confidence limits, and \eqn{p}-values.}
#'   \item{\code{residuals}}{the expected shortfall pseudo-response residuals,
#'   defined as the pseudo-response constructed using the final fitted
#'   conditional quantile values minus the final fitted expected shortfall
#'   values.}
#'   \item{\code{tau}}{the quantile level used for estimation.}
#'   \item{\code{method}}{the covariance estimation method, always
#'   \code{"rwb"} for random weighting bootstrap.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#'   \item{\code{level}}{the confidence level used for the reported confidence
#'   intervals.}
#'   \item{\code{B}}{the number of bootstrap replications.}
#'   \item{\code{cov}}{the estimated bootstrap covariance matrix, if
#'   \code{covariance = TRUE}.}
#' }
#'
#' @examples
#' set.seed(21)
#' dat <- data.frame(x = runif(1000))
#' dat$y <- 2 + dat$x + rnorm(1000)
#' fit <- irock(y ~ x, data = dat, tau = 0.1, bins = 3, ds = 0.02)
#' summary(fit, B = 200)
#' @export
summary.irock <- function(object, method = "rwb", covariance = FALSE,
                          level = 0.95, B = 500L, ...) {
  if (!inherits(object, "irock")) {
    stop("'object' must be an irock object.", call. = FALSE)
  }
  if (!identical(method, "rwb")) {
    stop("'method' must be \"rwb\".", call. = FALSE)
  }
  .irock_flag(covariance, "covariance")
  .irock_number(level, "level", upper = 1, closed = FALSE)
  .irock_count(B, "B")
  if (B <= 1L) stop("'B' must be greater than 1.", call. = FALSE)

  control <- attr(object, ".irock_control", exact = TRUE)
  if (is.null(control)) {
    stop("Refit the model with the updated irock() before requesting a summary.",
         call. = FALSE)
  }
  mf <- object$model
  x <- stats::model.matrix(object$terms, mf, contrasts.arg = object$contrasts)
  y <- stats::model.response(mf)
  theta_hat <- object$coefficients
  if (!identical(colnames(x), names(theta_hat))) {
    stop("Could not reconstruct the original model matrix.", call. = FALSE)
  }
  n <- nrow(x)
  parts <- .irock_predictors(mf[-attr(object$terms, "response")],
                             control$discrete, n)
  groups <- if (control$binning == "disjoint") {
    .irock_disjoint(parts, control$bins)
  } else {
    .irock_knn(x, control$neighbors, control$centers, control$medoids,
               representative = control$representative)
  }
  lower_tail <- object$tail == "lower"
  sign_y <- if (lower_tail) -1 else 1
  alpha <- if (lower_tail) 1 - object$tau else object$tau
  y <- sign_y * y

  coef.rwb <- matrix(NA_real_, B, length(theta_hat),
                     dimnames = list(NULL, names(theta_hat)))
  for (b in seq_len(B)) {
    w <- stats::rexp(n, rate = 1)
    w <- w / mean(w)
    fit.b <- tryCatch(
      .irock_fit(x, y, alpha, parts, control$delta, control$ds,
                 control$wins, control$bins, control$quantile_method,
                 control$ll_weight, control$binning, control$neighbors,
                 control$centers, control$medoids, groups = groups,
                 obs_weights = w, return_quantile = FALSE),
      error = function(e) {
        stop(sprintf("Random weighting bootstrap replication %d failed: %s",
                     b, conditionMessage(e)), call. = FALSE)
      }
    )
    coef.rwb[b, ] <- sign_y * fit.b$coefficients
  }

  centered <- sweep(coef.rwb, 2L, theta_hat, "-")
  cov_hat <- stats::cov(centered)
  se <- sqrt(diag(cov_hat))
  probs <- c((1 - level) / 2, (1 + level) / 2)
  limits <- apply(centered, 2L, stats::quantile, probs = probs, names = FALSE)
  p <- rep(NA_real_, length(se))
  positive <- se > 0
  p[positive] <- 2 * stats::pnorm(abs(theta_hat[positive] / se[positive]),
                                  lower.tail = FALSE)
  coef <- cbind(
    Estimate = theta_hat,
    "Std. Error" = se,
    "Lower CI" = theta_hat - limits[2L, ],
    "Upper CI" = theta_hat - limits[1L, ],
    "P-value" = p
  )
  ans <- list(call = object$call, terms = object$terms,
              coefficients.results = coef, residuals = object$residuals,
              tau = object$tau, method = method, tail = object$tail,
              level = level, B = B)
  if (covariance) ans$cov <- cov_hat
  class(ans) <- "summary.irock"
  ans
}
