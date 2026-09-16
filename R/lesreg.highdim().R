#' Fitting High-Dimensional Linear Expected Shortfall Regression Models
#'
#' @description
#' \code{lesreg.highdim} is used to fit high-dimensional linear expected
#' shortfall regression models by a two-step penalized procedure. It also
#' computes the associated penalized quantile regression fit used in the first
#' step of estimation.
#'
#' @usage
#' lesreg.highdim(x, y, alpha = 0.05, standardize = TRUE,
#'                tuning.method = c("CV", "BIC"), upperTail = FALSE)
#'
#' @param x numeric design matrix of dimension \eqn{n * p}, without an
#'   intercept column.
#' @param y vector of observations of length \eqn{n}, or a one-column matrix.
#' @param alpha the tail probability to be used. It must be a single number
#'   strictly between 0 and 1. This is the same tail level as \code{tau} in
#'   \code{\link{lesreg}}.
#' @param standardize logical. If \code{TRUE}, the covariates are standardized
#'   before fitting the penalized model sequence. The coefficients are returned
#'   on the original scale by the underlying fitting functions.
#' @param tuning.method character string specifying the method used to select
#'   the tuning parameters. Available choices are \code{"CV"} and \code{"BIC"}.
#' @param upperTail logical. If \code{FALSE}, the lower-tail expected shortfall
#'   regression is fitted. If \code{TRUE}, the upper-tail expected shortfall
#'   regression is fitted by applying the lower-tail procedure to the negative
#'   response and then flipping the signs of the resulting estimators.
#'
#' @details
#' This function implements a high-dimensional two-step estimator for linear
#' expected shortfall regression. By default, the lower-tail expected shortfall
#' regression model is fitted. The input matrix \code{x}
#' should not contain an intercept column. Intercepts for the quantile and
#' expected shortfall regressions are estimated separately.
#'
#' In the first step, a penalized quantile regression is fitted at level
#' \code{alpha}. If \code{tuning.method = "CV"}, the tuning parameter is selected
#' by cross-validation using \code{\link[conquer]{conquer.cv.reg}}. If
#' \code{tuning.method = "BIC"}, the tuning parameter is selected by a BIC-type
#' criterion using \code{\link[conquer]{conquer.reg}}.
#'
#' Let \eqn{\hat q} denote the fitted conditional quantile. The pseudo-response
#' \deqn{
#' z = \hat q + \alpha^{-1} 1(y \leq \hat q)(y - \hat q)
#' }
#' is then regressed on the same design matrix by penalized least squares using
#' \code{\link[glmnet]{glmnet}}. The resulting coefficients estimate the
#' high-dimensional linear expected shortfall function.
#'
#' If \code{upperTail = TRUE}, the upper-tail model is fitted by applying the
#' same lower-tail procedure to \code{-y}. The fitted quantile and expected
#' shortfall coefficients are then multiplied by \code{-1} and returned on the
#' original response scale.
#' Inference for a single coefficient can be obtained by
#' \code{\link{lesreg.highdim.inference}}.
#'
#' @return
#' \code{lesreg.highdim} returns an object of class
#' \code{"lesreg.highdim"}.
#'
#' An object of class \code{"lesreg.highdim"} is a list containing at least the
#' following components:
#' \describe{
#'   \item{theta0_hat}{intercept from the second-step expected shortfall
#'   regression fit.}
#'   \item{theta_hat}{slope coefficients from the second-step expected
#'   shortfall regression fit.}
#'   \item{beta0_hat}{intercept from the first-step quantile regression fit.}
#'   \item{beta_hat}{slope coefficients from the first-step quantile regression
#'   fit.}
#'   \item{resid}{residuals from the second-step fit, equal to the
#'   pseudo-response minus the fitted expected shortfall values.}
#'   \item{lambda_beta}{selected tuning parameter from the first-step quantile
#'   regression fit.}
#'   \item{lambda_theta}{selected tuning parameter from the second-step expected
#'   shortfall regression fit.}
#'   \item{tau}{the tail probability used for estimation.}
#'   \item{upperTail}{logical flag indicating whether the upper tail was fitted.}
#'   \item{tail}{character string indicating whether the fitted model is for the
#'   \code{"lower"} or \code{"upper"} tail.}
#'   \item{standardize}{logical flag indicating whether covariates were
#'   standardized in the penalized fitting routines.}
#'   \item{tuning.method}{the tuning parameter selection method used.}
#'   \item{call}{the matched call.}
#'   \item{y}{the original response vector.}
#'   \item{y.fit}{the response vector used internally for fitting. This is
#'   \code{y} for lower-tail fitting and \code{-y} for upper-tail fitting.}
#'   \item{theta_hat_fit}{second-step slope coefficients on the internal fitting
#'   scale.}
#'   \item{theta0_hat_fit}{second-step intercept on the internal fitting scale.}
#'   \item{beta_hat_fit}{first-step slope coefficients on the internal fitting
#'   scale.}
#'   \item{beta0_hat_fit}{first-step intercept on the internal fitting scale.}
#'   \item{resid_fit}{second-step residuals on the internal fitting scale.}
#' }
#'
#' @seealso
#' \code{\link{lesreg}}, \code{\link{lesreg.highdim.inference}},
#' \code{\link[conquer]{conquer.cv.reg}},
#' \code{\link[conquer]{conquer.reg}},
#' \code{\link[glmnet]{glmnet}},
#' \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#' set.seed(2026)
#' n <- 100
#' p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - x[, 2] + rnorm(n)
#'
#' fit <- lesreg.highdim(x, y, alpha = 0.05)
#' fit
#'
#' fit.bic <- lesreg.highdim(x, y, alpha = 0.05, tuning.method = "BIC")
#' fit.bic
#'
#' fit.upper <- lesreg.highdim(x, y, alpha = 0.05, upperTail = TRUE)
#' fit.upper
#' @export
lesreg.highdim <- function(x, y, alpha = 0.05, standardize = TRUE, tuning.method = c("CV", "BIC"), upperTail = FALSE){
  if(!is.logical(upperTail) || length(upperTail) != 1L || is.na(upperTail)){
    stop("'upperTail' must be TRUE or FALSE")
  }
  call <- match.call()
  tuning.method <- match.arg(tuning.method)
  if(is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  p <- ncol(x)
  y_original <- y
  if(upperTail){
    y_fit <- -y
  } else{
    y_fit <- y
  }
  if(tuning.method == "BIC"){
    quant_model_bic <- .lesreg_hd_quantile_bic(
      x,
      y_fit,
      alpha,
      standardize = standardize
    )
    lambda_beta <- quant_model_bic$lambda_select
    model_select <- quant_model_bic$model_select
    beta_hat <- model_select$coeff[-1]
    beta0 <- model_select$coeff[1]
    Z_2step <- apply(
      cbind(y_fit, x),
      MARGIN = 1,
      FUN = .lesreg_hd_Z_beta,
      beta = beta_hat,
      alpha = alpha,
      beta_0 = beta0
    )
    ES_model <- .lesreg_hd_LS_bic(
      x,
      Z_2step,
      standardize = standardize
    )
    lambda_theta <- ES_model$lambda_select
    theta_model <- ES_model$model_select
  } else{
    quan_model <- conquer::conquer.cv.reg(
      x,
      y_fit,
      tau = alpha,
      h = max(0.05, sqrt(alpha * (1 - alpha)) * (log(p) / n)^(1 / 4))
    )
    lambda_beta <- quan_model$lambda
    if(is.null(lambda_beta)){
      lambda_beta <- quan_model$lambda.min
    }
    beta_hat <- quan_model$coeff.min[-1]
    beta0 <- quan_model$coeff.min[1]
    Z_2step <- beta0 + x %*% beta_hat +
      (y_fit - beta0 - x %*% beta_hat) *
      ifelse(y_fit < beta0 + x %*% beta_hat, 1, 0) / alpha
    cv_model <- glmnet::cv.glmnet(
      x,
      Z_2step,
      standardize = standardize
    )
    lambda_theta <- cv_model$lambda.min
    theta_model <- glmnet::glmnet(
      x,
      Z_2step,
      lambda = lambda_theta,
      standardize = standardize
    )
  }
  theta_hat <- theta_model$beta
  theta0_hat <- theta_model$a0
  resid <- Z_2step - stats::predict(theta_model, newx = x)
  theta_hat_fit <- theta_hat
  theta0_hat_fit <- theta0_hat
  beta_hat_fit <- beta_hat
  beta0_hat_fit <- beta0
  resid_fit <- resid
  if(upperTail){
    theta_hat <- -theta_hat
    theta0_hat <- -theta0_hat
    beta_hat <- -beta_hat
    beta0 <- -beta0
    resid <- -resid
  }
  ans <- list(
    theta0_hat = theta0_hat,
    theta_hat = theta_hat,
    beta0_hat = beta0,
    beta_hat = beta_hat,
    resid = resid,
    lambda_beta = lambda_beta,
    lambda_theta = lambda_theta,
    tau = alpha,
    upperTail = upperTail,
    tail = if(upperTail) "upper" else "lower",
    standardize = standardize,
    tuning.method = tuning.method,
    call = call,
    y = y_original,
    y.fit = y_fit,
    theta_hat_fit = theta_hat_fit,
    theta0_hat_fit = theta0_hat_fit,
    beta_hat_fit = beta_hat_fit,
    beta0_hat_fit = beta0_hat_fit,
    resid_fit = resid_fit
  )
  class(ans) <- "lesreg.highdim"
  ans
}
