#' Fitting Linear Expected Shortfall Regression Models
#'
#' @description
#' \code{lesreg} is used to fit linear expected shortfall regression models.
#' It also computes the associated linear quantile regression fit used in the
#' first step of estimation.
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
#' @param tau the tail probability to be used. This is a single number strictly
#'   between 0 and 1. If \code{upperTail = FALSE}, \code{tau} is the lower-tail
#'   quantile level. If \code{upperTail = TRUE}, \code{tau} is the upper-tail
#'   probability and the corresponding quantile level is \eqn{1-\tau}.
#' @param data an optional data frame, list or environment containing the
#'   variables in the model. If not found in \code{data}, the variables are
#'   taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in the fitting process.
#' @param weights an optional vector of observation weights to be used in the
#'   fitting process. The weights must be positive. If supplied, weighted
#'   quantile regression and weighted least squares are used in the two steps.
#' @param na.action a function which indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} setting of
#'   \code{\link{options}}.
#' @param method the method to be used for the first-step quantile regression.
#'   This is passed to \code{\link[quantreg]{rq.fit}} or
#'   \code{\link[quantreg]{rq.wfit}}. If \code{method = "model.frame"}, the
#'   model frame is returned.
#' @param model,x,y,qr logicals. If \code{TRUE}, the corresponding components
#'   of the fit, the model frame, the model matrix, the response, and the QR
#'   decomposition, are returned.
#' @param contrasts an optional list. See the \code{contrasts.arg} argument of
#'   \code{\link{model.matrix.default}}.
#' @param eps a non-negative numerical tolerance used when constructing the
#'   tail indicator.
#' @param lm.tol tolerance passed to \code{\link[stats]{lm.fit}} or
#'   \code{\link[stats]{lm.wfit}} in the second-step least-squares fit.
#' @param singular.ok logical. If \code{FALSE}, a singular fit in the
#'   second-step least-squares regression is an error.
#' @param upperTail logical. If \code{FALSE}, the lower-tail expected shortfall
#'   is estimated. If \code{TRUE}, the upper-tail expected shortfall is
#'   estimated.
#' @param ... additional arguments to be passed to the low-level fitting
#'   functions, and ultimately to \code{\link[quantreg]{rq.fit}} or
#'   \code{\link[quantreg]{rq.wfit}}.
#'
#' @details
#' Models for \code{lesreg} are specified symbolically. A typical model has the
#' form \code{response ~ terms}, where \code{response} is the numeric response
#' vector and \code{terms} is a series of terms which specifies a linear
#' predictor for the conditional quantile and expected shortfall functions.
#'
#' A formula has an implied intercept term. To remove this use either
#' \code{y ~ x - 1} or \code{y ~ 0 + x}. See \code{\link{formula}} for more
#' details of allowed formulae.
#'
#' If \code{upperTail = FALSE}, the fitted model is
#' \deqn{
#' Q_\tau(Y \mid X) = X\beta,
#' }
#' together with
#' \deqn{
#' ES^L_\tau(Y \mid X)
#' =
#' E(Y \mid X, Y \leq Q_\tau(Y \mid X))
#' =
#' X\theta.
#' }
#' If \code{upperTail = TRUE}, the fitted expected shortfall model is
#' \deqn{
#' ES^U_\tau(Y \mid X)
#' =
#' E(Y \mid X, Y \geq Q_{1-\tau}(Y \mid X))
#' =
#' X\theta.
#' }
#'
#' All of \code{weights}, \code{subset} and \code{na.action} are evaluated in
#' the same way as variables in \code{formula}, that is first in \code{data}
#' and then in the environment of \code{formula}.
#'
#' \code{lesreg} calls the lower level functions \code{\link{lesreg.fit}} for
#' unweighted fitting and \code{\link{lesreg.wfit}} for weighted fitting. For
#' programming only, these functions may be called directly.
#'
#' @section Method:
#' The function computes a two-step estimate. First, a linear quantile
#' regression is fitted at level \code{tau}. Let \eqn{\hat q} denote the fitted
#' conditional quantile. The pseudo-response
#' \deqn{
#' w = \hat q + \tau^{-1} 1(y \leq \hat q + \epsilon)(y-\hat q)
#' }
#' is then regressed on the same model matrix by least squares. The resulting
#' coefficients estimate the linear expected shortfall function.
#'
#' For upper-tail fitting, the lower-tail procedure is applied internally to
#' the transformed response and model matrix. Fitted values are returned on the
#' original response scale.
#'
#' @return
#' \code{lesreg} returns an object of class \code{"lesreg"}.
#'
#' The functions \code{\link{summary.lesreg}} and \code{\link{print.lesreg}}
#' are used to obtain and print summaries of the results. The generic accessor
#' functions \code{\link{coef}}, \code{\link{fitted}} and
#' \code{\link{residuals}} can be used to extract useful features of the
#' fitted object.
#'
#' An object of class \code{"lesreg"} is a list containing at least the
#' following components:
#' \describe{
#'   \item{coefficients}{a named vector of expected shortfall regression
#'   coefficients.}
#'   \item{beta.hat}{a named vector of coefficients from the first-step
#'   quantile regression.}
#'   \item{theta.hat}{a named vector of coefficients from the second-step
#'   expected shortfall regression.}
#'   \item{residuals}{the residuals from the second-step regression, that is
#'   the pseudo-response minus fitted expected shortfall values.}
#'   \item{fitted.values}{the fitted expected shortfall values.}
#'   \item{quantile.residuals}{the residuals from the first-step quantile
#'   regression.}
#'   \item{quantile.fitted}{the fitted conditional quantile values.}
#'   \item{es.fitted}{the fitted conditional expected shortfall values.}
#'   \item{pseudo.response}{the pseudo-response used in the second-step
#'   least-squares regression.}
#'   \item{pseudo.residuals}{the residuals from the second-step least-squares
#'   regression.}
#'   \item{response.residuals}{the response minus fitted expected shortfall
#'   values.}
#'   \item{indicator}{the tail indicator used to construct the
#'   pseudo-response.}
#'   \item{rank}{the numeric rank of the fitted second-step linear model.}
#'   \item{weights}{the specified weights, or \code{NULL} for unweighted fits.}
#'   \item{df.residual}{the residual degrees of freedom.}
#'   \item{call}{the matched call.}
#'   \item{terms}{the \code{\link{terms}} object used.}
#'   \item{contrasts}{only where relevant, the contrasts used.}
#'   \item{xlevels}{only where relevant, a record of the levels of the factors
#'   used in fitting.}
#'   \item{method}{the method used for the first-step quantile regression.}
#'   \item{tau}{the tail probability used for estimation.}
#'   \item{tail}{a character string indicating whether the lower or upper tail
#'   was fitted.}
#'   \item{eps}{the numerical tolerance used in constructing the tail
#'   indicator.}
#'   \item{lm.tol}{the tolerance used in the second-step least-squares fit.}
#'   \item{singular.ok}{the value of \code{singular.ok} used in fitting.}
#'   \item{qr}{if requested, the QR decomposition from the second-step
#'   least-squares fit.}
#'   \item{y}{if requested, the response used.}
#'   \item{x}{if requested, the model matrix used.}
#'   \item{model}{if requested, the model frame used.}
#'   \item{na.action}{where relevant, information returned by
#'   \code{\link{model.frame}} on the special handling of \code{NA}s.}
#'   \item{quantile.fit}{the raw first-step quantile regression fit.}
#'   \item{es.fit}{the raw second-step least-squares fit.}
#' }
#'
#' In addition, fitted objects contain \code{x.fit} and \code{y.fit}, the model
#' matrix and response actually used by the low-level fitting routine.
#'
#' @seealso
#' \code{\link{summary.lesreg}}, \code{\link{lesreg.fit}},
#' \code{\link{lesreg.wfit}}
#'
#' \code{\link[quantreg]{rq}}, \code{\link[quantreg]{rq.fit}} and
#' \code{\link[quantreg]{rq.wfit}} for quantile regression.
#'
#' \code{\link[stats]{lm}}, \code{\link[stats]{lm.fit}} and
#' \code{\link[stats]{lm.wfit}} for least-squares regression.
#'
#' The generic functions \code{\link{coef}}, \code{\link{fitted}} and
#' \code{\link{residuals}}.
#'
#' @examples
#' set.seed(2026)
#' n <- 200
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' y <- 1 + 2 * x1 - x2 + rnorm(n)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#' fit <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05)
#' fit
#' summary(fit)
#' fit.upper <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05, upperTail = TRUE)
#' fit.upper
#' w <- stats::rexp(n)
#' fit.w <- lesreg(y ~ x1 + x2, data = dat, tau = 0.05, weights = w)
#' fit.w
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
  fit$eps <- eps
  fit$lm.tol <- lm.tol
  fit$singular.ok <- singular.ok
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
