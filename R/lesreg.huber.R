#' Robust Low-Dimensional Linear Expected Shortfall Regression
#'
#' @description
#' Fits a linear expected shortfall (ES) regression model in the
#' low-dimensional setting with potentially heavy-tailed response data using
#' the robust two-step framework of He et al. (2023). The first step fits the
#' corresponding linear quantile regression model (Koenker and Bassett, 1978),
#' and the second step applies the adaptive Huber regression procedure of
#' Pan and Zhou (2022) to the resulting pseudo-response.
#'
#' @usage
#' lesreg.huber(
#'   formula, tau = 0.05, data, subset,
#'   na.action, method = "br",
#'   model = TRUE, x = FALSE, y = FALSE,
#'   contrasts = NULL, eps = 1e-10,
#'   huber.epsilon = 1e-4,
#'   huber.maxit = 500L,
#'   upperTail = FALSE, ...
#' )
#'
#' @param formula an object of class \code{"formula"} specifying the model to
#'   be fitted. The current adaptive Huber implementation requires an
#'   intercept.
#' @param tau the quantile level associated with the ES to be estimated. This
#'   must be a single number strictly between 0 and 1. If
#'   \code{upperTail = FALSE}, the function fits the lower-tail ES associated
#'   with quantile level \code{tau}. If \code{upperTail = TRUE}, the function
#'   fits the upper-tail ES associated with quantile level \code{1 - tau}.
#'   See \code{\link{lesreg}} for details.
#' @param data an optional data frame, list, or environment containing the
#'   variables in the model. If a variable is not found in \code{data}, it is
#'   taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in fitting.
#' @param na.action a function specifying the action to be taken when the data
#'   contain \code{NA}s. The default is determined by the \code{na.action}
#'   setting of \code{\link{options}}.
#' @param method the method used for the first-step quantile regression. This
#'   is passed to \code{\link[quantreg]{rq.fit}}. If
#'   \code{method = "model.frame"}, the model frame is returned without
#'   fitting the model.
#' @param model,x,y logicals indicating whether the model frame, model matrix,
#'   and response, respectively, should be returned as components of the
#'   fitted object.
#' @param contrasts an optional list. See the \code{contrasts.arg} argument of
#'   \code{\link{model.matrix.default}}.
#' @param eps a non-negative numerical tolerance used when constructing the
#'   tail indicator.
#' @param huber.epsilon positive convergence tolerance used by the adaptive
#'   Huber regression algorithm.
#' @param huber.maxit positive integer giving the maximum number of iterations
#'   for the adaptive Huber regression algorithm.
#' @param upperTail logical. If \code{FALSE}, the lower-tail ES is estimated.
#'   If \code{TRUE}, the upper-tail ES is estimated.
#' @param ... additional arguments passed to the first-step quantile regression
#'   fitting routine.
#'
#' @details
#' \code{lesreg.huber} implements a robust version of the two-step linear ES
#' regression procedure described in \code{\link{lesreg}}. The first-step
#' quantile regression and pseudo-response \eqn{\hat{W}_i} construction follow the same
#' principle as in \code{lesreg}, whereas the second-step least-squares
#' regression is replaced by adaptive Huber regression to provide robustness
#' against heavy-tailed response data.
#'
#' Specifically, the ES regression coefficients are estimated by minimizing
#' the empirical Huber loss:
#' \deqn{
#' \hat{\boldsymbol{\theta}} = \mathrm{argmin}_{\boldsymbol{\theta}} \frac{1}{N}
#' \sum_{i = 1}^nH_{r}(\hat{W}_i - \boldsymbol{X}_i^{\intercal}\boldsymbol{\theta}),
#' }
#' where \eqn{H_r(x) = x^2 \mathbf{1}(|x|\leq r)/2 + (r|x| - r^2/2)\mathbf{1}(|x|>r)}
#' denotes the Huber loss function (Huber, 1973), and the robustification
#' parameter \eqn{r} is adaptively selected using the procedure of Pan and
#' Zhou (2022).
#'
#' If \code{upperTail = TRUE}, the upper-tail model is fitted using the same
#' transformation principle as in \code{\link{lesreg}}. Fitted quantities are
#' returned on the original response scale. The current implementation
#' requires a model containing an intercept.
#'
#' @return
#' \code{lesreg.huber} returns an object of class \code{"lesreg.huber"}.
#' Important components include:
#' \describe{
#'   \item{\code{coefficients}}{a named vector of adaptive Huber expected
#'   shortfall regression coefficients.}
#'   \item{\code{beta.hat}}{a named vector of linear quantile regression
#'   coefficients from the first-step quantile regression.}
#'   \item{\code{fitted.values}}{the fitted expected shortfall values.}
#'   \item{\code{quantile.fitted}}{the fitted conditional quantile values.}
#'   \item{\code{pseudo.response}}{the pseudo-response used in the second-step
#'   adaptive Huber regression.}
#'   \item{\code{huber.r}}{the adaptively selected Huber robustification
#'   parameter.}
#'   \item{\code{huber.iteration}}{the number of adaptive Huber iterations
#'   performed.}
#'   \item{\code{tau}}{the tail probability used for estimation.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#'   \item{\code{call}}{the call that produced the fitted object.}
#' }
#'
#' @references
#' He, X., Tan, K. M., and Zhou, W.-X. (2023).
#' Robust estimation and inference for expected shortfall regression with many
#' regressors.
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology},
#' 85(4), 1223--1246.
#'
#' Huber, P. J. (1973).
#' Robust estimation: asymptotics, conjectures and Monte Carlo.
#' \emph{The Annals of Statistics}, 1, 799--821.
#'
#' Koenker, R., and Bassett, G. Jr. (1978).
#' Regression quantiles.
#' \emph{Econometrica}, 46(1), 33--50.
#'
#' Pan, X. and Zhou, W.-X. (2022).
#' \emph{adaHuber: Adaptive Huber estimation and regression}.
#' R package version 1.1.
#' Available at \url{https://github.com/XiaoouPan/adaHuber}.
#'
#' @seealso
#' \code{\link{lesreg.huber.fit}},
#' \code{\link{lesreg}},
#' \code{\link{lesreg.fit}}
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
#'     stats::integrate(stats::qt, df = df, lower = 0, upper = tau)$value / tau
#'
#' theta.true.lower <- c(1 + est, 2, -1)
#' theta.true.upper <- c(1 - est, 2, -1)
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit lower-quantile robust ES
#' fit <- lesreg.huber(
#'   y ~ x1 + x2,
#'   data = dat,
#'   tau = tau
#' )
#' tab.lower <- rbind(fit$coefficients, theta.true.lower)
#' rownames(tab.lower) <- c('theta.est.lower', 'theta.true.lower')
#' print(tab.lower)
#'
#' # Fit upper-quantile robust ES
#' fit.upper <- lesreg.huber(
#' y ~ x1 + x2,
#'   data = dat,
#'   tau = 0.1,
#'   upperTail = TRUE
#' )
#' tab.upper <- rbind(fit.upper$coefficients, theta.true.upper)
#' rownames(tab.upper) <- c('theta.est.upper', 'theta.true.upper')
#' print(tab.upper)
#' @export

lesreg.huber <- function(formula, tau = 0.05, data, subset,
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
  fit$pseudo.response <- drop(
    fit$pseudo.response
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
