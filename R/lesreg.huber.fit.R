#' Low-Level Fitter for Robust Linear Expected Shortfall Regression
#'
#' @description
#' Low-level fitting function for robust linear expected shortfall (ES)
#' regression using the adaptive Huber two-step procedure implemented by
#' \code{\link{lesreg.huber}}. The function works directly with a numeric
#' design matrix and response vector.
#'
#' For formula-based model specification, data handling, subsets,
#' missing-value handling, and factor variables, use
#' \code{\link{lesreg.huber}}.
#'
#' @usage
#' lesreg.huber.fit(
#'   x, y, tau = 0.05, method = "br",
#'   upperTail = FALSE, eps = 1e-10,
#'   huber.epsilon = 1e-4, huber.maxit = 500L, ...
#' )
#'
#' @param x numeric design matrix of dimension \eqn{n \times p}. The current
#'   implementation requires exactly one intercept column.
#' @param y numeric response vector of length \eqn{n}, or a one-column numeric
#'   matrix.
#' @param tau the quantile level associated with the ES to be estimated. This
#'   must be a single number strictly between 0 and 1. If
#'   \code{upperTail = FALSE}, the function fits the lower-tail ES associated
#'   with quantile level \code{tau}. If \code{upperTail = TRUE}, the function
#'   fits the upper-tail ES associated with quantile level \code{1 - tau}.
#'   See \code{\link{lesreg.huber}} for details.
#' @param method the method used for the first-step quantile regression. This
#'   is passed to \code{\link[quantreg]{rq.fit}}.
#' @param upperTail logical. If \code{FALSE}, the lower-tail ES is estimated.
#'   If \code{TRUE}, the upper-tail ES is estimated.
#' @param eps a non-negative numerical tolerance used when constructing the
#'   tail indicator.
#' @param huber.epsilon positive numerical convergence tolerance used by the
#'   adaptive Huber regression algorithm.
#' @param huber.maxit positive integer giving the maximum number of iterations
#'   for the adaptive Huber regression algorithm.
#' @param ... additional arguments passed to
#'   \code{\link[quantreg]{rq.fit}}.
#'
#' @details
#' \code{lesreg.huber.fit} is the low-level fitting function used by
#' \code{\link{lesreg.huber}} for unweighted robust linear ES regression. It
#' does not perform formula processing or construct a model frame and
#' therefore requires the design matrix and response to be supplied directly.
#'
#' The function uses the same first-step quantile regression and
#' pseudo-response construction as \code{\link{lesreg.fit}}, but replaces the
#' second-step least-squares regression with adaptive Huber regression. The
#' Huber robustification parameter is selected internally and returned as
#' \code{huber.r}.
#'
#' If \code{upperTail = TRUE}, the upper-tail model is fitted using the same
#' transformation principle as in \code{\link{lesreg.huber}}. Fitted
#' quantities are returned on the original response scale, whereas
#' \code{x.fit} and \code{y.fit} retain the internal fitting scale.
#'
#' The current implementation requires exactly one intercept column in
#' \code{x}.
#'
#' @return
#' A list containing at least the following components:
#' \describe{
#'   \item{\code{coefficients}}{a named vector of adaptive Huber expected
#'   shortfall regression coefficients.}
#'   \item{\code{beta.hat}}{a named vector of linear quantile regression
#'   coefficients from the first-step quantile regression.}
#'   \item{\code{x.fit}}{the design matrix used internally for fitting.}
#'   \item{\code{y.fit}}{the response used internally for fitting.}
#'   \item{\code{residuals}}{the residuals from the second-step adaptive Huber
#'   regression, defined as the pseudo-response minus the fitted expected
#'   shortfall values.}
#'   \item{\code{fitted.values}}{the fitted expected shortfall values.}
#'   \item{\code{quantile.residuals}}{the residuals from the first-step
#'   quantile regression.}
#'   \item{\code{quantile.fitted}}{the fitted conditional quantile values.}
#'   \item{\code{pseudo.response}}{the pseudo-response used in the second-step
#'   adaptive Huber regression.}
#'   \item{\code{response.residuals}}{the observed response minus the fitted
#'   expected shortfall values.}
#'   \item{\code{indicator}}{the tail indicator used to construct the
#'   pseudo-response.}
#'   \item{\code{rank}}{the numeric rank of the design matrix used in
#'   fitting.}
#'   \item{\code{df.residual}}{the residual degrees of freedom based on the
#'   rank of the design matrix.}
#'   \item{\code{huber.r}}{the adaptively selected Huber robustification
#'   parameter.}
#'   \item{\code{huber.iteration}}{the number of adaptive Huber iterations
#'   performed.}
#'   \item{\code{huber.epsilon}}{the Huber convergence tolerance used.}
#'   \item{\code{huber.maxit}}{the maximum number of Huber iterations allowed.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#' }
#'
#' @seealso
#' \code{\link{lesreg.huber}},
#' \code{\link{lesreg.fit}},
#' \code{\link{lesreg}}
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
#' est <-
#'     stats::integrate(stats::qt, df = df, lower = 0, upper = tau)$value/tau
#' theta.true.lower <- c(1 + est, 2, -1)
#'
#' x <- model.matrix(~ x1 + x2)
#'
#' fit <- lesreg.huber.fit(
#'   x,
#'   y,
#'   tau = 0.1
#' )
#' print(theta.true.lower)
#' print(fit$coefficients)
#' print(fit$huber.r)
#'
#' @export

lesreg.huber.fit <- function(x, y, tau = 0.05, method = "br",
                             upperTail = FALSE, eps = 1e-10,
                             huber.epsilon = 1e-4,
                             huber.maxit = 500L, ...){
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
  if(!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0 || tau >= 1){
    stop("'tau' must be a single number strictly between 0 and 1")
  }
  if(!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps < 0){
    stop("'eps' must be a non-negative numeric scalar")
  }
  if(!is.numeric(huber.epsilon) || length(huber.epsilon) != 1L || is.na(huber.epsilon) || huber.epsilon <= 0){
    stop("'huber.epsilon' must be a positive numeric scalar")
  }
  if(length(huber.maxit) != 1L || is.na(huber.maxit) || huber.maxit < 1){
    stop("'huber.maxit' must be a positive integer")
  }
  huber.maxit <- as.integer(huber.maxit)
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
      rank = 0L,
      df.residual = length(y),
      assign = attr(x, "assign"),
      qr = NULL,
      quantile.fit = NULL,
      es.fit = NULL,
      huber.r = NA_real_,
      huber.iteration = 0L,
      tail = if (upperTail) "upper" else "lower",
      contrasts = attr(x, "contrasts")
    ))
  }
  x.fit <- x
  y.fit <- y
  if(upperTail){
    x.fit <- -x.fit
    y.fit <- -y.fit
  }
  fit.qr <- quantreg::rq.fit(
    x.fit,
    y.fit,
    tau = tau,
    method = method,
    ...
  )
  beta.hat <- fit.qr$coefficients
  q.hat.fit <- drop(x.fit %*% beta.hat)
  quantile.residuals.fit <- drop(y.fit - q.hat.fit)
  indicator <- as.numeric(y.fit <= q.hat.fit + eps)
  w.fit <- drop(q.hat.fit + (1 / tau) * indicator * (y.fit - q.hat.fit))
  assign.x <- attr(x, "assign")
  intercept.col <- integer(0L)
  if(!is.null(assign.x)){
    intercept.col <- which(assign.x == 0L)
  }
  if(length(intercept.col) == 0L && !is.null(colnames(x)) && "(Intercept)" %in% colnames(x)){
    intercept.col <- match("(Intercept)", colnames(x))
  }
  if(length(intercept.col) != 1L){
    stop("The current adaptive Huber implementation requires exactly one intercept column.")
  }
  huber.x <- x.fit[, -intercept.col, drop = FALSE]
  fit.es <- adaHuberReg(
    X = huber.x,
    Y = w.fit,
    epsilon = huber.epsilon,
    iteMax = huber.maxit
  )
  theta.hat <- numeric(p)
  non.intercept <- setdiff(seq_len(p), intercept.col)
  theta.hat[intercept.col] <- fit.es$coef[1L] / x.fit[1L, intercept.col]
  theta.hat[non.intercept] <- fit.es$coef[-1L]
  es.hat.fit <- drop(x.fit %*% theta.hat)
  q.hat <- if(upperTail){
    -q.hat.fit
  } else{
    q.hat.fit
  }
  w <- if(upperTail){
    -w.fit
  } else{
    w.fit
  }
  es.hat <- if(upperTail){
    -es.hat.fit
  } else{
    es.hat.fit
  }
  quantile.residuals <- y - q.hat
  pseudo.residuals <- w - es.hat
  response.residuals <- y - es.hat
  dn <- colnames(x)
  if(is.null(dn)){
    dn <- paste0("x", seq_len(p))
  }
  names(beta.hat) <- dn
  names(theta.hat) <- dn
  qr.x <- qr(x.fit)
  rank.x <- qr.x$rank
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
    rank = rank.x,
    df.residual = n - rank.x,
    assign = attr(x, "assign"),
    qr = qr.x,
    quantile.fit = fit.qr,
    es.fit = fit.es,
    huber.r = fit.es$tau,
    huber.iteration = fit.es$iteration,
    huber.epsilon = huber.epsilon,
    huber.maxit = huber.maxit,
    tail = if(upperTail) "upper" else "lower"
  )
  fit$contrasts <- attr(x, "contrasts")
  fit
}
