#' Fitter Function for Linear Expected Shortfall Regression
#'
#' @description
#' Basic computing engine used by \code{\link{lesreg}} to fit an unweighted
#' linear expected shortfall regression model.
#'
#' This function works directly with a numeric design matrix and response
#' vector. It does not process formulas, model frames, missing values,
#' subsets, weights, or factor levels. These tasks are handled by the
#' higher-level \code{\link{lesreg}} function.
#'
#' @usage
#' lesreg.fit(x, y, tau = 0.05, method = "br",
#'            upperTail = FALSE, eps = 1e-10,
#'            lm.tol = 1e-7, singular.ok = TRUE, ...)
#'
#' @param x Design matrix of dimension \eqn{n * p}.
#' @param y Vector of observations of length \eqn{n}.
#' @param tau Tail probability used in estimation. It must be a single number
#'   strictly between 0 and 1. The default is \code{0.05}. If
#'   \code{upperTail = FALSE}, \code{tau} is the lower-tail quantile level. If
#'   \code{upperTail = TRUE}, \code{tau} represents the upper-tail probability,
#'   so the corresponding quantile level is \eqn{1-\tau}.
#' @param method Method used by \code{\link[quantreg]{rq.fit}} in the
#'   first-step quantile regression. The default is \code{"br"}.
#' @param upperTail Logical. If \code{FALSE}, the lower-tail expected shortfall
#'   is estimated. If \code{TRUE}, the upper-tail expected shortfall is
#'   estimated by applying the lower-tail procedure to the transformed data
#'   \eqn{(-y, -x)}. The default is \code{FALSE}.
#' @param eps Non-negative numerical tolerance used when constructing the tail
#'   indicator. The default is \code{1e-10}.
#' @param lm.tol Tolerance passed to the \code{tol} argument of
#'   \code{\link[stats]{lm.fit}} in the second-step least-squares regression.
#'   The default is \code{1e-7}.
#' @param singular.ok Logical. If \code{FALSE}, a singular fit in the
#'   second-step least-squares regression is an error.
#' @param ... Additional arguments passed to \code{\link[quantreg]{rq.fit}}.
#'
#' @details
#' \code{lesreg.fit} is the low-level unweighted fitting function for the
#' two-step estimator of linear expected shortfall regression.
#'
#' If \code{upperTail = FALSE}, the function estimates the lower-tail model
#' \deqn{
#' ES^L_\tau(Y \mid X)
#' =
#' E{Y \mid X, Y \leq Q_\tau(Y \mid X)}
#' =
#' X\theta.
#' }
#'
#' In the first step, a linear conditional quantile regression model is fitted
#' by \code{\link[quantreg]{rq.fit}}:
#' \deqn{
#' Q_\tau(Y \mid X)=X\beta.
#' }
#'
#' Let \eqn{\hat{q}=X\hat{\beta}} denote the fitted conditional quantile. The
#' pseudo-response is then constructed as
#' \deqn{
#' w = \hat{q} + \frac{1}{\tau}1(y \leq \hat{q})(y-\hat{q}).
#' }
#'
#' In the second step, \code{\link[stats]{lm.fit}} is used to regress the
#' pseudo-response \eqn{w} on the same design matrix \eqn{X}. The resulting
#' coefficient vector estimates \eqn{\theta}.
#'
#' If \code{upperTail = TRUE}, the function estimates the upper-tail model
#' \deqn{
#' ES^U_\tau(Y \mid X)
#' =
#' E{Y \mid X, Y \geq Q_{1-\tau}(Y \mid X)}
#' =
#' X\theta.
#' }
#'
#' This is computed by defining
#' \deqn{
#' Y^*=-Y,\qquad X^*=-X.
#' }
#' Since
#' \deqn{
#' Y \geq Q_{1-\tau}(Y \mid X)
#' \iff
#' Y^* \leq Q_\tau(Y^* \mid X^*),
#' }
#' the upper-tail problem for \eqn{(Y,X)} is transformed into a lower-tail
#' problem for \eqn{(Y^*,X^*)}. The same lower-tail procedure is then applied
#' to the transformed data. The estimated coefficient vector is reported as
#' \eqn{\theta}, and fitted values are converted back to the original response
#' scale.
#'
#' This function is intended for advanced use. Most users should call
#' \code{\link{lesreg}} instead.
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{coefficients}{Estimated expected shortfall regression coefficients.
#'   Same as \code{theta.hat}.}
#'   \item{beta.hat}{Estimated coefficients from the first-step conditional
#'   quantile regression.}
#'   \item{theta.hat}{Estimated coefficients from the second-step conditional
#'   expected shortfall regression.}
#'   \item{residuals}{Residuals from the second-step regression of the
#'   pseudo-response on \code{x}. Same as \code{pseudo.residuals}.}
#'   \item{fitted.values}{Fitted conditional expected shortfall values on the
#'   original response scale. Same as \code{es.fitted}.}
#'   \item{quantile.residuals}{Residuals from the first-step quantile
#'   regression, computed on the original response scale.}
#'   \item{quantile.fitted}{Fitted conditional quantile values on the original
#'   response scale. For \code{upperTail = TRUE}, these correspond to the
#'   fitted \eqn{1-\tau} conditional quantile.}
#'   \item{es.fitted}{Fitted conditional expected shortfall values on the
#'   original response scale.}
#'   \item{pseudo.response}{Pseudo-response used in the second-step
#'   least-squares regression, reported on the original response scale.}
#'   \item{pseudo.residuals}{Residuals from the second-step regression of the
#'   pseudo-response on \code{x}.}
#'   \item{response.residuals}{Difference between the original response
#'   \code{y} and the fitted expected shortfall values.}
#'   \item{indicator}{Tail indicator vector used to construct the
#'   pseudo-response. For \code{upperTail = FALSE}, this corresponds to the
#'   lower-tail event. For \code{upperTail = TRUE}, this corresponds to the
#'   upper-tail event on the original response scale.}
#'   \item{rank}{Integer giving the rank of the second-step design matrix.}
#'   \item{df.residual}{Degrees of freedom of residuals from the second-step
#'   regression.}
#'   \item{assign}{Assignment attribute from the design matrix.}
#'   \item{qr}{QR decomposition from the second-step least-squares regression.}
#'   \item{quantile.fit}{Raw fit object returned by the first-step
#'   \code{\link[quantreg]{rq.fit}} call.}
#'   \item{es.fit}{Raw fit object returned by the second-step
#'   \code{\link[stats]{lm.fit}} call.}
#'   \item{tail}{Character string indicating the fitted tail, either
#'   \code{"lower"} or \code{"upper"}.}
#'   \item{contrasts}{Contrasts attribute from the design matrix, if present.}
#' }
#'
#' @seealso
#' \code{\link{lesreg}}, \code{\link{lesreg.wfit}},
#' \code{\link[quantreg]{rq.fit}}, \code{\link[stats]{lm.fit}}
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' x <- model.matrix(~ x1 + x2)
#'
#' ## Lower-tail fit
#' fit.lower <- lesreg.fit(x, y, tau = 0.05)
#' fit.lower$beta.hat
#' fit.lower$theta.hat
#' head(fit.lower$quantile.fitted)
#' head(fit.lower$es.fitted)
#'
#' ## Upper-tail fit
#' fit.upper <- lesreg.fit(x, y, tau = 0.05, upperTail = TRUE)
#' fit.upper$beta.hat
#' fit.upper$theta.hat
#' head(fit.upper$quantile.fitted)
#' head(fit.upper$es.fitted)
#'
#' @export

lesreg.fit <- function(x, y, tau = 0.05, method = "br", upperTail = FALSE, eps = 1e-10, lm.tol = 1e-7, singular.ok = TRUE, ...){
  if(is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  if(n == 0L)
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if(!is.logical(upperTail) || length(upperTail) != 1L || is.na(upperTail))
    stop("'upperTail' must be TRUE or FALSE")
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
    if (upperTail) {
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
  fit.qr <- quantreg::rq.fit(x.fit, y.fit, tau = tau, method = method, ...)
  beta.hat <- fit.qr$coefficients
  quantile.residuals.fit <- drop(fit.qr$residuals)
  q.hat.fit <- drop(y.fit - quantile.residuals.fit)
  indicator <- as.numeric(y.fit <= q.hat.fit + eps)
  w.fit <- drop(q.hat.fit + (1 / tau) * indicator * (y.fit - q.hat.fit))
  fit.es <- stats::lm.fit(x = x.fit, y = w.fit, tol = lm.tol, singular.ok = singular.ok)
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

