#' Fitter Function for Adaptive Huber Expected Shortfall Regression
#'
#' @description
#' This is the basic computing engine for unweighted linear expected shortfall
#' regression with an adaptive Huber second-step estimator. It works directly
#' with a numeric design matrix and response vector.
#'
#' For formula processing, missing-value handling, subsets, and factor
#' variables, use \code{\link{lesreg.huber}}.
#'
#' @usage
#' lesreg.huber.fit(
#'   x, y, tau = 0.05, method = "br",
#'   upperTail = FALSE, eps = 1e-10,
#'   huber.epsilon = 1e-4, huber.maxit = 500L, ...
#' )
#'
#' @param x numeric design matrix of dimension \eqn{n * p}. The current
#'   implementation requires exactly one intercept column.
#' @param y numeric response vector of length \eqn{n}, or a one-column matrix.
#' @param tau tail probability. It must be a single number strictly between
#'   0 and 1.
#' @param method method used for the first-step quantile regression; passed to
#'   \code{\link[quantreg]{rq.fit}}.
#' @param upperTail logical. If \code{FALSE}, the lower-tail expected shortfall
#'   is fitted. If \code{TRUE}, the upper-tail expected shortfall is fitted.
#' @param eps non-negative numerical tolerance used when constructing the tail
#'   indicator.
#' @param huber.epsilon positive numerical convergence tolerance used by the
#'   adaptive Huber regression algorithm.
#' @param huber.maxit positive integer giving the maximum number of adaptive
#'   Huber iterations.
#' @param ... additional arguments passed to
#'   \code{\link[quantreg]{rq.fit}}.
#'
#' @details
#' The estimator is computed in two steps. First, linear quantile regression
#' is fitted at level \code{tau}. Let \eqn{\hat q_i} denote the fitted
#' conditional quantile. The pseudo-response is
#'
#' \deqn{
#' w_i =
#' \hat q_i +
#' \tau^{-1}
#' 1(y_i \leq \hat q_i + \epsilon)
#' (y_i-\hat q_i).
#' }
#'
#' Unlike \code{\link{lesreg.fit}}, which fits the second step by least
#' squares, \code{lesreg.huber.fit} regresses the pseudo-response on the
#' covariates using adaptive Huber regression.
#'
#' The adaptive Huber tuning parameter is estimated internally and is returned
#' as \code{huber.tau}.
#'
#' For upper-tail estimation, the corresponding lower-tail problem is fitted
#' internally after a sign transformation, and fitted quantities are returned
#' on the original response scale.
#'
#' The current implementation requires an intercept because the underlying
#' adaptive Huber fitting routine constructs its own intercept internally.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{coefficients}{adaptive Huber expected shortfall regression
#'     coefficients.}
#'   \item{beta.hat}{first-step quantile regression coefficients.}
#'   \item{theta.hat}{second-step adaptive Huber expected shortfall
#'     coefficients.}
#'   \item{x.fit}{design matrix used in fitting.}
#'   \item{y.fit}{response used internally in fitting.}
#'   \item{residuals}{second-step pseudo-response residuals.}
#'   \item{fitted.values}{fitted expected shortfall values.}
#'   \item{quantile.residuals}{first-step quantile residuals.}
#'   \item{quantile.fitted}{fitted conditional quantile values.}
#'   \item{es.fitted}{fitted conditional expected shortfall values.}
#'   \item{pseudo.response}{pseudo-response used in the adaptive Huber
#'     second step.}
#'   \item{pseudo.residuals}{pseudo-response residuals.}
#'   \item{response.residuals}{response minus fitted expected shortfall values.}
#'   \item{indicator}{tail indicator used to construct the pseudo-response.}
#'   \item{rank}{rank of the design matrix.}
#'   \item{df.residual}{residual degrees of freedom based on the design matrix.}
#'   \item{quantile.fit}{raw first-step quantile regression fit.}
#'   \item{es.fit}{raw adaptive Huber fit returned by the internal C++ engine.}
#'   \item{huber.tau}{adaptively selected Huber robustification parameter.}
#'   \item{huber.iteration}{number of adaptive Huber iterations performed.}
#'   \item{huber.epsilon}{Huber convergence tolerance used.}
#'   \item{huber.maxit}{maximum number of Huber iterations allowed.}
#'   \item{tail}{either \code{"lower"} or \code{"upper"}.}
#' }
#'
#' @seealso
#' \code{\link{lesreg.huber}},
#' \code{\link{lesreg.fit}},
#' \code{\link{lesreg}}
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#'
#' x <- model.matrix(~ x1 + x2)
#'
#' fit <- lesreg.huber.fit(
#'   x,
#'   y,
#'   tau = 0.1
#' )
#'
#' fit$theta.hat
#' fit$huber.tau
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
      huber.tau = NA_real_,
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
    rank = rank.x,
    df.residual = n - rank.x,
    assign = attr(x, "assign"),
    qr = qr.x,
    quantile.fit = fit.qr,
    es.fit = fit.es,
    huber.tau = fit.es$tau,
    huber.iteration = fit.es$iteration,
    huber.epsilon = huber.epsilon,
    huber.maxit = huber.maxit,
    tail = if(upperTail) "upper" else "lower"
  )
  fit$contrasts <- attr(x, "contrasts")
  fit
}
