#' Low-Level Weighted Fitter for Linear Expected Shortfall Regression
#'
#' @description
#' Low-level fitting function for weighted linear expected shortfall (ES)
#' regression using the two-step regression approach implemented by
#' \code{\link{lesreg}}. The function works directly with a numeric design
#' matrix, response vector, and observation weights.
#'
#' For formula-based model specification, data handling, subsets,
#' missing-value handling, and factor variables, use \code{\link{lesreg}}.
#'
#' @usage
#' lesreg.wfit(x, y, weights = NULL, tau = 0.05, method = "br",
#'             upperTail = FALSE, eps = 1e-10,
#'             lm.tol = 1e-7, singular.ok = TRUE, ...)
#'
#' @param x numeric design matrix of dimension \eqn{n \times p}.
#' @param y numeric response vector of length \eqn{n}, or a one-column numeric
#'   matrix.
#' @param weights a vector of positive observation weights of length \eqn{n}
#'   to be used in fitting. If \code{NULL}, unit weights are used.
#' @param tau the quantile level associated with the ES to be estimated. This
#'   must be a single number strictly between 0 and 1. If
#'   \code{upperTail = FALSE}, the function fits the lower-tail ES associated
#'   with quantile level \code{tau}. If \code{upperTail = TRUE}, the function
#'   fits the upper-tail ES associated with quantile level \code{1 - tau}.
#'   See \code{\link{lesreg}} for details.
#' @param method the method used for the first-step quantile regression. This
#'   is passed to \code{\link[quantreg]{rq.wfit}}.
#' @param upperTail logical. If \code{FALSE}, the lower-tail ES is estimated.
#'   If \code{TRUE}, the upper-tail ES is estimated.
#' @param eps a non-negative numerical tolerance used when constructing the
#'   tail indicator.
#' @param lm.tol numerical tolerance passed to \code{\link[stats]{lm.wfit}} in
#'   the second-step weighted least-squares fit.
#' @param singular.ok logical. If \code{FALSE}, a singular fit in the
#'   second-step weighted least-squares regression results in an error.
#' @param ... additional arguments passed to
#'   \code{\link[quantreg]{rq.wfit}}.
#'
#' @details
#' \code{lesreg.wfit} is the low-level weighted fitting function used by
#' \code{\link{lesreg}}. It does not perform formula processing or construct a
#' model frame and therefore requires the design matrix and response to be
#' supplied directly.
#'
#' The function implements the weighted version of the two-step estimation
#' procedure described in \code{\link{lesreg}} and \code{\link{lesreg.fit}}.
#' In the first step, a weighted linear quantile regression is fitted using
#' \code{\link[quantreg]{rq.wfit}}. In the second step, the pseudo-response
#' constructed from the fitted conditional quantile is regressed on the same
#' design matrix using weighted least squares through
#' \code{\link[stats]{lm.wfit}}.
#'
#' If \code{upperTail = TRUE}, the lower-tail procedure is applied internally
#' to the transformed response and design matrix. Fitted values are returned
#' on the original response scale. See \code{\link{lesreg}} for the
#' statistical formulation of the lower- and upper-tail ES regression models.
#'
#' If \code{x} has no columns, the function returns empty coefficient vectors
#' without calling the first-step or second-step fitting routines.
#'
#' @return
#' A list containing at least the following components:
#' \describe{
#'   \item{\code{coefficients}}{a named vector of expected shortfall regression
#'   coefficients.}
#'   \item{\code{beta.hat}}{a named vector of linear quantile regression
#'   coefficients from the first-step quantile regression.}
#'   \item{\code{x.fit}}{the design matrix used by the low-level fitting
#'   routine.}
#'   \item{\code{y.fit}}{the response used by the low-level fitting routine.}
#'   \item{\code{residuals}}{the residuals from the second-step weighted
#'   least-squares regression, defined as the pseudo-response minus the fitted
#'   expected shortfall values.}
#'   \item{\code{fitted.values}}{the fitted expected shortfall values.}
#'   \item{\code{quantile.residuals}}{the residuals from the first-step quantile
#'   regression.}
#'   \item{\code{quantile.fitted}}{the fitted conditional quantile values.}
#'   \item{\code{pseudo.response}}{the pseudo-response used in the second-step
#'   weighted least-squares regression.}
#'   \item{\code{response.residuals}}{the observed response minus the fitted
#'   expected shortfall values.}
#'   \item{\code{indicator}}{the vector of ``hitting indicators" that compute
#'   \eqn{\mathbf{1}(Y_i \leq \boldsymbol{X}_i^{\intercal}\hat{\boldsymbol{\beta}})}
#'   when \code{upperTail = FALSE}, and
#'   \eqn{\mathbf{1}(Y_i \geq \boldsymbol{X}_i^{\intercal}\hat{\boldsymbol{\beta}})}
#'   when \code{upperTail = TRUE}.}
#'   \item{\code{weights}}{the vector of weights used in fitting.}
#'   \item{\code{rank}}{the numeric rank of the fitted second-step linear model.}
#'   \item{\code{df.residual}}{the residual degrees of freedom.}
#'   \item{\code{assign}}{the assignment attribute of the design matrix.}
#'   \item{\code{qr}}{the QR decomposition from the second-step weighted
#'   least-squares fit.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#'   \item{\code{contrasts}}{the contrasts attribute of the design matrix, if
#'   present.}
#' }
#'
#' @seealso
#' \code{\link{lesreg}} for the user-level interface and a detailed
#' description of the two-step ES regression procedure;
#' \code{\link{lesreg.fit}} for the corresponding unweighted low-level fitter.
#'
#' \code{\link[quantreg]{rq.wfit}} for the first-step weighted quantile
#' regression and \code{\link[stats]{lm.wfit}} for the second-step weighted
#' least-squares regression.
#' @examples
#' set.seed(2026)
#' n <- 500
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' tau <- 0.05
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' x <- model.matrix(~ x1 + x2)
#' weights <- stats::rexp(n)
#' theta.true.lower <- c(1 +
#'     stats::integrate(stats::qnorm, lower = 0, upper = tau)$value/tau, 2, -1)
#' theta.true.upper <- c(1 +
#'     stats::integrate(stats::qnorm, lower = 1-tau, upper = 1)$value/tau, 2, -1 )
#'
#' ## Fit weighted linear ES regression for lower tail
#' fit <- lesreg.wfit(x, y, weights = weights, tau = tau)
#' fit$coefficients
#' theta.true.lower
#'
#' ## Fit weighted linear ES regression for upper tail
#' fit.upper <- lesreg.wfit(x, y, weights = weights, tau = tau, upperTail = TRUE)
#' fit.upper$coefficients
#' theta.true.upper
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
      x.fit = x.fit,
      y.fit = y.fit,
      residuals = y,
      fitted.values = 0 * y,
      quantile.residuals = y,
      quantile.fitted = 0 * y,
      pseudo.response = y,
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
    residuals = pseudo.residuals,
    fitted.values = es.hat,
    quantile.residuals = quantile.residuals,
    quantile.fitted = q.hat,
    pseudo.response = w,
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
