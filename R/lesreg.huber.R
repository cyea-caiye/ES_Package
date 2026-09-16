#' Fitting Adaptive Huber Expected Shortfall Regression Models
#'
#' @description
#' \code{lesreg.huber} fits linear expected shortfall regression models using
#' quantile regression in the first step and adaptive Huber regression in the
#' second step.
#'
#' @usage
#' lesreg.huber(
#'   formula, tau = 0.05, data, subset, weights,
#'   na.action, method = "br",
#'   model = TRUE, x = FALSE, y = FALSE,
#'   contrasts = NULL, eps = 1e-10,
#'   huber.epsilon = 1e-4,
#'   huber.maxit = 500L,
#'   upperTail = FALSE, ...
#' )
#'
#' @param formula an object of class \code{"formula"} specifying the response
#'   and linear predictor.
#' @param tau tail probability. It must be a single number strictly between
#'   0 and 1.
#' @param data optional data frame, list, or environment containing the
#'   variables in the model.
#' @param subset optional vector specifying observations to be used.
#' @param weights observation weights. Weighted adaptive Huber expected
#'   shortfall regression is not currently implemented, and supplying this
#'   argument produces an error.
#' @param na.action function specifying the handling of missing values.
#' @param method method used for the first-step quantile regression; passed to
#'   \code{\link[quantreg]{rq.fit}}. If \code{method = "model.frame"},
#'   the model frame is returned.
#' @param model logical. If \code{TRUE}, the model frame is retained.
#' @param x logical. If \code{TRUE}, the model matrix is retained.
#' @param y logical. If \code{TRUE}, the response is retained.
#' @param contrasts optional list of contrasts passed to
#'   \code{\link[stats]{model.matrix}}.
#' @param eps non-negative numerical tolerance used when constructing the tail
#'   indicator.
#' @param huber.epsilon positive convergence tolerance used by the adaptive
#'   Huber regression algorithm.
#' @param huber.maxit positive integer giving the maximum number of adaptive
#'   Huber iterations.
#' @param upperTail logical. If \code{FALSE}, lower-tail expected shortfall is
#'   estimated. If \code{TRUE}, upper-tail expected shortfall is estimated.
#' @param ... additional arguments passed to the first-step quantile regression
#'   fitting routine.
#'
#' @details
#' A linear quantile regression is first fitted at the specified tail
#' probability. Let \eqn{\hat q_i} be the resulting fitted conditional
#' quantile. The pseudo-response
#'
#' \deqn{
#' w_i =
#' \hat q_i +
#' \tau^{-1}
#' 1(y_i \leq \hat q_i + \epsilon)
#' (y_i-\hat q_i)
#' }
#'
#' is then fitted using adaptive Huber regression.
#'
#' Thus, \code{lesreg.huber} uses the same first-step quantile regression and
#' pseudo-response construction as \code{\link{lesreg}}, but replaces the
#' least-squares second step with adaptive Huber regression.
#'
#' For upper-tail fitting, an equivalent lower-tail problem is fitted
#' internally after a sign transformation.
#'
#' The current implementation requires a model containing an intercept.
#' Weighted adaptive Huber expected shortfall regression is not currently
#' implemented.
#'
#' @return
#' An object of class \code{"lesreg.huber"}. Important components include:
#' \describe{
#'   \item{coefficients}{adaptive Huber expected shortfall coefficients.}
#'   \item{beta.hat}{first-step quantile regression coefficients.}
#'   \item{theta.hat}{adaptive Huber expected shortfall coefficients.}
#'   \item{fitted.values}{fitted expected shortfall values.}
#'   \item{quantile.fitted}{fitted conditional quantile values.}
#'   \item{pseudo.response}{pseudo-response used in the second step.}
#'   \item{huber.tau}{adaptively selected Huber robustification parameter.}
#'   \item{huber.iteration}{number of Huber iterations performed.}
#'   \item{tau}{tail probability used in estimation.}
#'   \item{tail}{either \code{"lower"} or \code{"upper"}.}
#'   \item{call}{matched function call.}
#' }
#'
#' @seealso
#' \code{\link{lesreg.huber.fit}},
#' \code{\link{lesreg}},
#' \code{\link{lesreg.fit}}
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#'
#' # Heavy-tailed errors
#' eps <- rt(n, df = 2.5)
#' y <- 1 + 2 * x1 - x2 + eps
#'
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' fit <- lesreg.huber(
#'   y ~ x1 + x2,
#'   data = dat,
#'   tau = 0.1
#' )
#' fit
#'
#' fit.upper <- lesreg.huber(
#'   y ~ x1 + x2,
#'   data = dat,
#'   tau = 0.1,
#'   upperTail = TRUE
#' )
#' fit.upper
#'
#' @export

lesreg.huber <- function(formula, tau = 0.05, data, subset, weights,
                         na.action, method = "br",
                         model = TRUE, x = FALSE, y = FALSE,
                         contrasts = NULL, eps = 1e-10,
                         huber.epsilon = 1e-4,
                         huber.maxit = 500L,
                         upperTail = FALSE, ...) {
  if(!inherits(formula, "formula")){
    stop("'formula' must be a formula")
  }
  if(!is.logical(model) || length(model) != 1L || is.na(model)){
    stop("'model' must be TRUE or FALSE")
  }
  if(!is.logical(x) || length(x) != 1L || is.na(x)){
    stop("'x' must be TRUE or FALSE")
  }
  if(!is.logical(y) || length(y) != 1L || is.na(y)){
    stop("'y' must be TRUE or FALSE")
  }
  if(!is.character(method) || length(method) != 1L || is.na(method)){
    stop("'method' must be a single character string")
  }
  if(!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0 || tau >= 1){
    stop("'tau' must be a single number strictly between 0 and 1")
  }
  if(!is.logical(upperTail) || length(upperTail) != 1L || is.na(upperTail)){
    stop("'upperTail' must be TRUE or FALSE")
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
  ret.x <- x
  ret.y <- y
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(
    c("formula", "data", "subset", "na.action"),
    names(mf),
    0L
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(method == "model.frame"){
    return(mf)
  }
  mt <- attr(mf, "terms")
  if(attr(mt, "intercept") != 1L){
    stop(
      "the current adaptive Huber implementation requires a model with an intercept"
    )
  }
  Y <- stats::model.response(mf)
  if(NCOL(Y) != 1L){
    stop("'y' must be a numeric vector or a one-column matrix")
  }

  if(is.matrix(Y)){
    Y <- drop(Y)
  }

  if(!is.numeric(Y)){
    stop("'y' must be numeric")
  }
  X <- stats::model.matrix(
    mt,
    mf,
    contrasts
  )
  fit <- lesreg.huber.fit(
    x = X,
    y = Y,
    tau = tau,
    method = method,
    upperTail = upperTail,
    eps = eps,
    huber.epsilon = huber.epsilon,
    huber.maxit = huber.maxit,
    ...
  )
  class(fit) <- "lesreg.huber"
  fit$na.action <- attr(mf, "na.action")
  fit$formula <- formula
  fit$terms <- mt
  fit$xlevels <- stats::.getXlevels(mt, mf)
  fit$call <- call
  fit$tau <- tau
  fit$weights <- NULL
  fit$method <- method
  fit$eps <- eps
  fit$huber.epsilon <- huber.epsilon
  fit$huber.maxit <- huber.maxit
  fit$residuals <- drop(
    fit$residuals
  )
  fit$fitted.values <- drop(
    fit$fitted.values
  )
  fit$quantile.residuals <- drop(
    fit$quantile.residuals
  )
  fit$quantile.fitted <- drop(
    fit$quantile.fitted
  )
  fit$es.fitted <- drop(
    fit$es.fitted
  )
  fit$pseudo.response <- drop(
    fit$pseudo.response
  )
  fit$pseudo.residuals <- drop(
    fit$pseudo.residuals
  )
  fit$response.residuals <- drop(
    fit$response.residuals
  )
  if(model){
    fit$model <- mf
  }
  if(ret.x){
    fit$x <- X
  }
  if(ret.y){
    fit$y <- Y
  }
  fit
}
