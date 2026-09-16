#' Weighted Fitter Function for Linear Expected Shortfall Regression
#'
#' @description
#' This is the basic computing engine called by \code{\link{lesreg}} for
#' weighted linear expected shortfall regression. It works directly with a
#' numeric design matrix, response vector, and observation weights.
#'
#' For formula processing, missing-value handling, subsets, and factor
#' variables, use \code{\link{lesreg}}.
#'
#' @usage
#' lesreg.wfit(x, y, weights = NULL, tau = 0.05, method = "br",
#'             upperTail = FALSE, eps = 1e-10,
#'             lm.tol = 1e-7, singular.ok = TRUE, ...)
#'
#' @param x design matrix of dimension \eqn{n * p}.
#' @param y vector of observations of length \eqn{n}, or a one-column matrix.
#' @param weights vector of positive weights of length \eqn{n} to be used in
#'   the fitting process. If \code{NULL}, unit weights are used.
#' @param tau the tail probability to be used. It must be a single number
#'   strictly between 0 and 1.
#' @param method method of computation for the first-step weighted quantile
#'   regression; passed to \code{\link[quantreg]{rq.wfit}}.
#' @param upperTail logical. If \code{FALSE}, the lower tail is fitted. If
#'   \code{TRUE}, the upper tail is fitted.
#' @param eps non-negative numerical tolerance used in constructing the tail
#'   indicator.
#' @param lm.tol tolerance passed to \code{\link[stats]{lm.wfit}} in the
#'   second-step weighted least-squares fit.
#' @param singular.ok logical. If \code{FALSE}, a singular second-step fit is
#'   an error.
#' @param ... optional arguments passed to \code{\link[quantreg]{rq.wfit}}.
#'
#' @details
#' The function computes a two-step weighted fit. First, a weighted linear
#' quantile regression is fitted by \code{\link[quantreg]{rq.wfit}}. The fitted
#' quantile is then used to form the pseudo-response
#' \deqn{
#' w = \hat q + \tau^{-1} 1(y \leq \hat q + \epsilon)(y-\hat q),
#' }
#' which is regressed on the same design matrix by
#' \code{\link[stats]{lm.wfit}} using the same weights.
#'
#' If \code{upperTail = TRUE}, the same lower-tail computation is applied to
#' the transformed data \eqn{-y} and \eqn{-x}. Fitted values are returned on the
#' original response scale.
#'
#' Fits without any columns return empty coefficient vectors and do not call
#' the first-step or second-step fitting routines.
#'
#' @return
#' a list with components
#' \describe{
#'   \item{coefficients}{\eqn{p} vector of expected shortfall regression coefficients.}
#'   \item{beta.hat}{\eqn{p} vector of first-step weighted quantile regression coefficients.}
#'   \item{theta.hat}{\eqn{p} vector of second-step weighted expected shortfall
#'   regression coefficients.}
#'   \item{x.fit}{design matrix used in fitting.}
#'   \item{y.fit}{response used in fitting.}
#'   \item{residuals}{\eqn{n} vector of second-step residuals.}
#'   \item{fitted.values}{\eqn{n} vector of fitted expected shortfall values.}
#'   \item{quantile.residuals}{\eqn{n} vector of first-step quantile residuals.}
#'   \item{quantile.fitted}{\eqn{n} vector of fitted quantile values.}
#'   \item{es.fitted}{\eqn{n} vector of fitted expected shortfall values.}
#'   \item{pseudo.response}{\eqn{n} vector used as the response in the
#'   second-step weighted least-squares fit.}
#'   \item{pseudo.residuals}{\eqn{n} vector of residuals from the second-step fit.}
#'   \item{response.residuals}{\eqn{n} vector of response residuals.}
#'   \item{indicator}{\eqn{n} vector indicating the observations used in the tail.}
#'   \item{weights}{\eqn{n} vector of weights used in fitting.}
#'   \item{rank}{integer, giving the rank of the second-step design matrix.}
#'   \item{df.residual}{degrees of freedom of residuals.}
#'   \item{assign}{the assignment attribute of the design matrix.}
#'   \item{qr}{the QR decomposition from the second-step fit.}
#'   \item{quantile.fit}{the object returned by the first-step
#'   \code{\link[quantreg]{rq.wfit}} call.}
#'   \item{es.fit}{the object returned by the second-step
#'   \code{\link[stats]{lm.wfit}} call.}
#'   \item{tail}{character string, either \code{"lower"} or \code{"upper"}.}
#'   \item{contrasts}{the contrasts attribute of the design matrix, if any.}
#' }
#'
#' @seealso
#' \code{\link{lesreg}}, \code{\link{lesreg.fit}},
#' \code{\link[quantreg]{rq.wfit}}, \code{\link[stats]{lm.wfit}}
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' x <- model.matrix(~ x1 + x2)
#' weights <- stats::rexp(n)
#' fit <- lesreg.wfit(x, y, weights = weights, tau = 0.05)
#' fit$theta.hat
#' fit.unit <- lesreg.wfit(x, y, tau = 0.05)
#' fit.unit$weights
#' fit.upper <- lesreg.wfit(x, y, weights = weights, tau = 0.05, upperTail = TRUE)
#' fit.upper$theta.hat
#'
#' @export


lesreg.wfit <- function(x, y, weights = NULL, tau = 0.05, method = "br",
                        upperTail = FALSE, eps = 1e-10, lm.tol = 1e-7,
                        singular.ok = TRUE, ...){
  if(is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  if(n == 0L)
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if(!is.logical(upperTail) || length(upperTail) != 1L || is.na(upperTail)){
    stop("'upperTail' must be TRUE or FALSE")
  }
  ny <- NCOL(y)
  if(ny != 1L){
    stop("'y' must be a numeric vector or a one-column matrix")
  }
  if(is.matrix(y)){
    y <- drop(y)
  }
  if(NROW(y) != n){
    stop("incompatible dimensions")
  }
  if(!is.numeric(y)){
    stop("'y' must be numeric")
  }
  if(is.null(weights)){
    weights <- rep(1, n)
  } else{
    weights <- as.vector(weights)
    if(length(weights) != n){
      stop("'weights' must have the same length as 'y'")
    }
    if(!is.numeric(weights)){
      stop("'weights' must be numeric")
    }
    if(any(!is.finite(weights))){
      stop("'weights' must be finite")
    }
    if(any(weights <= 0)){
      stop("'weights' must be positive")
    }
  }
  if(!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0 || tau >= 1){
    stop("'tau' must be a single number strictly between 0 and 1")
  }
  if(!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0){
    stop("'eps' must be a non-negative numeric scalar")
  }
  if(!is.numeric(lm.tol) || length(lm.tol) != 1L || is.na(lm.tol) || lm.tol < 0){
    stop("'lm.tol' must be a non-negative numeric scalar")
  }
  if(p == 0L){
    x.fit <- x
    y.fit <- y
    if(upperTail){
      x.fit <- -x.fit
      y.fit <- -y.fit
    }
    return(list(
      coefficients = numeric(),
      beta.hat = numeric(),
      theta.hat = numeric(),
      x.fit = x.fit,
      y.fit = y.fit,
      residuals = y,
      fitted.values = 0 * y,
      quantile.residuals = y,
      quantile.fitted = 0 * y,
      es.fitted = 0 * y,
      pseudo.response = y,
      pseudo.residuals = y,
      response.residuals = y,
      indicator = numeric(length(y)),
      weights = weights,
      rank = 0L,
      df.residual = length(y),
      assign = attr(x, "assign"),
      qr = NULL,
      quantile.fit = NULL,
      es.fit = NULL,
      tail = if(upperTail) "upper" else "lower",
      contrasts = attr(x, "contrasts")
    ))
  }
  x.fit <- x
  y.fit <- y
  if(upperTail){
    x.fit <- -x.fit
    y.fit <- -y.fit
  }
  fit.qr <- quantreg::rq.wfit(x = x.fit, y = y.fit, tau = tau, weights = weights, method = method, ...)
  beta.hat <- fit.qr$coefficients
  quantile.residuals.fit <- drop(fit.qr$residuals)
  q.hat.fit <- drop(y.fit - quantile.residuals.fit)
  indicator <- as.numeric(y.fit <= q.hat.fit + eps)
  w.fit <- drop(q.hat.fit + (1 / tau) * indicator * (y.fit - q.hat.fit))
  fit.es <- stats::lm.wfit(x = x.fit, y = w.fit, w = weights, tol = lm.tol, singular.ok = singular.ok)
  theta.hat <- fit.es$coefficients
  es.hat.fit <- fit.es$fitted.values
  q.hat <- if(upperTail) -q.hat.fit else q.hat.fit
  w <- if(upperTail) -w.fit else w.fit
  es.hat <- if(upperTail) -es.hat.fit else es.hat.fit
  quantile.residuals <- y - q.hat
  pseudo.residuals <- w - es.hat
  response.residuals <- y - es.hat
  dn <- colnames(x)
  if(is.null(dn)){
    dn <- paste0("x", seq_len(p))
  }
  names(beta.hat) <- dn
  names(theta.hat) <- dn
  fit <- list(
    coefficients = theta.hat,
    x.fit = x.fit,
    y.fit = y.fit,
    beta.hat = beta.hat,
    theta.hat = theta.hat,
    residuals = pseudo.residuals,
    fitted.values = es.hat,
    quantile.residuals = quantile.residuals,
    quantile.fitted = q.hat,
    es.fitted = es.hat,
    pseudo.response = w,
    pseudo.residuals = pseudo.residuals,
    response.residuals = response.residuals,
    indicator = indicator,
    weights = weights,
    rank = fit.es$rank,
    df.residual = fit.es$df.residual,
    assign = attr(x, "assign"),
    qr = fit.es$qr,
    quantile.fit = fit.qr,
    es.fit = fit.es,
    tail = if(upperTail) "upper" else "lower"
  )
  fit$contrasts <- attr(x, "contrasts")
  fit
}
