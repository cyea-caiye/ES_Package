#' Fitting Linear Expected Shortfall Regression Models
#'
#' @description
#' \code{lesreg} is used to fit linear conditional expected shortfall
#' regression models by a two-step procedure.
#'
#' It can be used to estimate a linear conditional quantile function and the
#' corresponding linear conditional expected shortfall function at a specified
#' tail probability \eqn{\tau}.
#'
#' @usage
#' lesreg(formula, tau = 0.05, data, subset, weights, na.action,
#'        method = "br", model = TRUE, x = FALSE, y = FALSE,
#'        qr = TRUE, contrasts = NULL, eps = 1e-10,
#'        lm.tol = 1e-7, singular.ok = TRUE,
#'        upperTail = FALSE, ...)
#'
#' @param formula an object of class \code{"formula"}: a symbolic description
#'   of the model to be fitted. The details of model specification are given
#'   under \code{Details}.
#' @param tau Tail probability used in estimation. It must be a single numeric
#'   value strictly between 0 and 1. The default is \code{0.05}. If
#'   \code{upperTail = FALSE}, \code{tau} is the lower-tail quantile level. If
#'   \code{upperTail = TRUE}, \code{tau} represents the upper-tail probability,
#'   so the corresponding quantile level is \eqn{1-\tau}.
#' @param data an optional data frame, list or environment containing the
#'   variables in the model. If not found in \code{data}, the variables are
#'   taken from the environment of \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights an optional vector of positive observation weights to be used
#'   in the fitting process. If supplied, \code{lesreg} calls
#'   \code{\link{lesreg.wfit}}. If omitted, \code{lesreg} calls
#'   \code{\link{lesreg.fit}} for an unweighted fit.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \code{\link{options}}.
#' @param method the method to be used for the first-step quantile regression.
#'   This argument is passed to \code{\link[quantreg]{rq.fit}} for unweighted
#'   fits and to \code{\link[quantreg]{rq.wfit}} for weighted fits. The default
#'   is \code{"br"}. If \code{method = "model.frame"}, the model frame is
#'   returned and no model is fitted.
#' @param model,x,y,qr logicals. If \code{TRUE}, the corresponding components
#'   of the fit are returned: the model frame, the model matrix, the response,
#'   and the QR decomposition from the second-step least-squares regression.
#' @param contrasts an optional list. See the \code{contrasts.arg} argument of
#'   \code{\link{model.matrix.default}}.
#' @param eps a non-negative numeric tolerance used when constructing the tail
#'   indicator. The default is \code{1e-10}.
#' @param lm.tol tolerance passed to the \code{tol} argument of the second-step
#'   least-squares regression. The default is \code{1e-7}.
#' @param singular.ok logical. If \code{FALSE}, a singular fit in the
#'   second-step least-squares regression is an error.
#' @param upperTail logical. If \code{FALSE}, the lower-tail expected shortfall
#'   is estimated. If \code{TRUE}, the upper-tail expected shortfall is
#'   estimated by applying the lower-tail procedure to the transformed data
#'   \eqn{(-Y,-X)}. The default is \code{FALSE}.
#' @param ... additional arguments to be passed to the low-level fitting
#'   function, and ultimately to \code{\link[quantreg]{rq.fit}} or
#'   \code{\link[quantreg]{rq.wfit}}.
#'
#' @details
#' Models for \code{lesreg} are specified symbolically. A typical model has the
#' form \code{response ~ terms}, where \code{response} is the numeric response
#' vector and \code{terms} is a series of terms which specifies a linear
#' predictor for the conditional quantile and expected shortfall functions.
#'
#' A formula has an implied intercept term. To remove this, use either
#' \code{y ~ x - 1} or \code{y ~ 0 + x}. See \code{\link{formula}} for more
#' details of allowed formulae.
#'
#' If \code{upperTail = FALSE}, the fitted model is based on
#' \deqn{
#' Q_\tau(Y \mid X) = X\beta
#' }
#' and
#' \deqn{
#' ES^L_\tau(Y \mid X)
#' =
#' E{Y \mid X, Y \leq Q_\tau(Y \mid X)}
#' =
#' X\theta.
#' }
#'
#' The estimation proceeds in two steps. First, \code{lesreg} constructs the
#' model frame and model matrix, and calls either \code{\link{lesreg.fit}} or
#' \code{\link{lesreg.wfit}} to fit the first-step linear quantile regression.
#' Let
#' \deqn{
#' \hat{q}=X\hat{\beta}
#' }
#' denote the fitted conditional quantile.
#'
#' Second, the pseudo-response
#' \deqn{
#' w = \hat{q} + \frac{1}{\tau}1(y \leq \hat{q})(y-\hat{q})
#' }
#' is constructed. The expected shortfall coefficient vector \eqn{\theta} is
#' then estimated by regressing \eqn{w} on the same design matrix \eqn{X}.
#'
#' If \code{weights} is not supplied, \code{lesreg} calls
#' \code{\link{lesreg.fit}} and performs unweighted estimation. If
#' \code{weights} is supplied, \code{lesreg} calls \code{\link{lesreg.wfit}}
#' and uses the supplied observation weights in both estimation steps.
#'
#' If \code{upperTail = TRUE}, \code{lesreg} estimates
#' \deqn{
#' ES^U_\tau(Y \mid X)
#' =
#' E{Y \mid X, Y \geq Q_{1-\tau}(Y \mid X)}
#' =
#' X\theta.
#' }
#' This case is computed internally by applying the lower-tail procedure to the
#' transformed response and design matrix
#' \deqn{
#' Y^*=-Y,\qquad X^*=-X.
#' }
#' The fitted coefficients are reported for the original upper-tail model, and
#' fitted values are reported on the original response scale.
#'
#' \code{lesreg} calls the lower-level function \code{\link{lesreg.fit}} for
#' unweighted fits and \code{\link{lesreg.wfit}} for weighted fits. For
#' programming with numeric model matrices directly, these lower-level fitting
#' functions may be used instead.
#'
#' All of \code{weights}, \code{subset} and \code{na.action} are evaluated in
#' the same way as variables in \code{formula}, that is first in \code{data}
#' and then in the environment of \code{formula}.
#'
#' @return
#' \code{lesreg} returns an object of class \code{"lesreg"}.
#'
#' An object of class \code{"lesreg"} is a list containing at least the
#' following components:
#' \describe{
#'   \item{coefficients}{a named vector of expected shortfall regression
#'   coefficients. Same as \code{theta.hat}.}
#'   \item{beta.hat}{a named vector of coefficients from the first-step
#'   conditional quantile regression.}
#'   \item{theta.hat}{a named vector of coefficients from the second-step
#'   conditional expected shortfall regression.}
#'   \item{residuals}{the residuals from the second-step regression of the
#'   pseudo-response on the model matrix. Same as \code{pseudo.residuals}.}
#'   \item{fitted.values}{the fitted conditional expected shortfall values on
#'   the original response scale. Same as \code{es.fitted}.}
#'   \item{quantile.residuals}{the residuals from the first-step quantile
#'   regression, computed on the original response scale.}
#'   \item{quantile.fitted}{the fitted conditional quantile values on the
#'   original response scale. For \code{upperTail = TRUE}, these correspond to
#'   the fitted \eqn{1-\tau} conditional quantile.}
#'   \item{es.fitted}{the fitted conditional expected shortfall values on the
#'   original response scale.}
#'   \item{pseudo.response}{the pseudo-response used in the second-step
#'   least-squares regression, reported on the original response scale.}
#'   \item{pseudo.residuals}{the residuals from the second-step regression of
#'   the pseudo-response on the model matrix.}
#'   \item{response.residuals}{the original response minus the fitted expected
#'   shortfall values.}
#'   \item{indicator}{the tail indicator used to construct the pseudo-response.
#'   For \code{upperTail = FALSE}, this corresponds to the lower-tail event. For
#'   \code{upperTail = TRUE}, this corresponds to the upper-tail event on the
#'   original response scale.}
#'   \item{rank}{the numeric rank of the second-step fitted linear model.}
#'   \item{weights}{the specified weights, or \code{NULL} for unweighted fits.}
#'   \item{df.residual}{the residual degrees of freedom from the second-step
#'   regression.}
#'   \item{call}{the matched call.}
#'   \item{terms}{the \code{\link{terms}} object used.}
#'   \item{contrasts}{only where relevant, the contrasts used.}
#'   \item{xlevels}{only where relevant, a record of the levels of the factors
#'   used in fitting.}
#'   \item{method}{the method used for the first-step quantile regression.}
#'   \item{tau}{the tail probability used for estimation.}
#'   \item{tail}{character string indicating the fitted tail, either
#'   \code{"lower"} or \code{"upper"}.}
#'   \item{qr}{if requested, the QR decomposition from the second-step
#'   least-squares regression.}
#'   \item{y}{if requested, the response used.}
#'   \item{x}{if requested, the model matrix used.}
#'   \item{model}{if requested, the model frame used.}
#'   \item{na.action}{where relevant, information returned by
#'   \code{\link{model.frame}} on the special handling of \code{NA}s.}
#'   \item{quantile.fit}{the raw fit object returned by the first-step
#'   quantile regression call.}
#'   \item{es.fit}{the raw fit object returned by the second-step
#'   least-squares regression call.}
#' }
#'
#' @seealso
#' \code{\link{lesreg.fit}} and \code{\link{lesreg.wfit}} for the underlying
#' low-level fitting functions.
#'
#' \code{\link[quantreg]{rq}}, \code{\link[quantreg]{rq.fit}}, and
#' \code{\link[quantreg]{rq.wfit}} for quantile regression.
#'
#' \code{\link[stats]{lm}}, \code{\link[stats]{lm.fit}}, and
#' \code{\link[stats]{lm.wfit}} for linear least-squares regression.
#'
#' The generic functions \code{\link{coef}}, \code{\link{residuals}},
#' \code{\link{fitted}}, and \code{\link{predict}}.
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Lower-tail unweighted fit
#' fit.lower <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05)
#' fit.lower$beta.hat
#' fit.lower$theta.hat
#' head(fit.lower$quantile.fitted)
#' head(fit.lower$es.fitted)
#'
#' ## Upper-tail unweighted fit
#' fit.upper <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05,
#'                     upperTail = TRUE)
#' fit.upper$tail
#' head(fit.upper$quantile.fitted)
#' head(fit.upper$es.fitted)
#'
#' ## Weighted fit
#' weights <- stats::rexp(n, rate = 1)
#' fit.w <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05,
#'                 weights = weights)
#' fit.w$beta.hat
#' fit.w$theta.hat
#'
#' ## Return model matrix and response
#' fit2 <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05,
#'                x = TRUE, y = TRUE)
#' names(fit2)
#'
#' ## Inspect the model frame without fitting
#' mf <- lesreg(y ~ x1 + x2, data = dat, method = "model.frame")
#' head(mf)
#'
#' @export

lesreg <- function(formula, tau = 0.05, data, subset, weights, na.action,
                   method = "br", model = TRUE, x = FALSE, y = FALSE,
                   qr = TRUE, contrasts = NULL, eps = 1e-10, lm.tol = 1e-7,
                   singular.ok = TRUE, upperTail = FALSE, ...){
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
  if(!is.logical(qr) || length(qr) != 1L || is.na(qr)){
    stop("'qr' must be TRUE or FALSE")
  }
  if(!is.character(method) || length(method) != 1L || is.na(method)){
    stop("'method' must be a single character string")
  }
  ret.x <- x
  ret.y <- y
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if(method == "model.frame"){
    return(mf)
  }
  mt <- attr(mf, "terms")
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
  weights <- as.vector(stats::model.weights(mf))
  if(!is.null(weights)){
    if(!is.numeric(weights)){
      stop("'weights' must be a numeric vector")
    }
    if(length(weights) != NROW(Y)){
      stop("'weights' must have the same length as the response")
    }
    if(any(!is.finite(weights))){
      stop("'weights' must be finite")
    }
    if(any(weights <= 0)){
      stop("'weights' must be positive")
    }
  }
  if(length(tau) != 1L || !is.numeric(tau) || is.na(tau) || tau <= 0 || tau >= 1){
    stop("'tau' must be a single number strictly between 0 and 1")
  }
  if(!is.logical(upperTail) || length(upperTail) != 1L || is.na(upperTail)){
    stop("'upperTail' must be TRUE or FALSE")
  }
  X <- stats::model.matrix(mt, mf, contrasts)
  if(is.null(weights)){
    fit <- lesreg.fit(
      x = X,
      y = Y,
      tau = tau,
      method = method,
      upperTail = upperTail,
      eps = eps,
      lm.tol = lm.tol,
      singular.ok = singular.ok,
      ...
    )
  } else{
    fit <- lesreg.wfit(
      x = X,
      y = Y,
      weights = weights,
      tau = tau,
      method = method,
      upperTail = upperTail,
      eps = eps,
      lm.tol = lm.tol,
      singular.ok = singular.ok,
      ...
    )
  }
  class(fit) <- "lesreg"
  fit$na.action <- attr(mf, "na.action")
  fit$formula <- formula
  fit$terms <- mt
  fit$xlevels <- stats::.getXlevels(mt, mf)
  fit$call <- call
  fit$tau <- tau
  fit$weights <- weights
  fit$method <- method
  fit$residuals <- drop(fit$residuals)
  fit$fitted.values <- drop(fit$fitted.values)
  fit$quantile.residuals <- drop(fit$quantile.residuals)
  fit$quantile.fitted <- drop(fit$quantile.fitted)
  fit$es.fitted <- drop(fit$es.fitted)
  fit$pseudo.response <- drop(fit$pseudo.response)
  fit$pseudo.residuals <- drop(fit$pseudo.residuals)
  fit$response.residuals <- drop(fit$response.residuals)
  if(model){
    fit$model <- mf
  }
  if(ret.x) fit$x <- X
  if(ret.y) fit$y <- Y
  if(!qr){
    fit$qr <- NULL
    if(!is.null(fit$es.fit)){
      fit$es.fit$qr <- NULL
    }
  }
  fit
}
