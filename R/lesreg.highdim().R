#' Fitting High-Dimensional Linear Expected Shortfall Regression Models
#'
#' @description
#' Fits high-dimensional linear expected shortfall (ES) regression models using
#' the two-step penalized estimation procedure of Zhang et al. (2025). The
#' function also fits the associated penalized quantile regression model used
#' in the first step of estimation.
#'
#' @param x numeric design matrix of dimension \eqn{n \times p}, without an
#'   intercept column.
#' @param y numeric response vector of length \eqn{n}, or a one-column numeric
#'   matrix.
#' @param alpha the quantile level associated with the ES to be estimated. This
#'   must be a single number strictly between 0 and 1. If
#'   \code{upperTail = FALSE}, the function fits the lower-tail ES associated
#'   with quantile level \code{alpha}. If \code{upperTail = TRUE}, the function
#'   fits the upper-tail ES associated with quantile level \code{1 - alpha}.
#'   This corresponds to \code{tau} in \code{\link{lesreg}}.
#' @param standardize logical. If \code{TRUE}, standardization is used in the
#'   penalized fitting procedures where supported. The coefficients are
#'   returned on the original covariate scale by the underlying fitting
#'   functions.
#' @param tuning.method character string specifying the method used to select
#'   the tuning parameters in the two estimation steps. Available choices are
#'   \code{"CV"} for cross-validation and \code{"BIC"} for BIC-type tuning
#'   parameter selection.
#' @param upperTail logical. If \code{FALSE}, the lower-tail ES is estimated.
#'   If \code{TRUE}, the upper-tail ES is estimated by applying the lower-tail
#'   procedure to the negative response and transforming the resulting
#'   estimators back to the original response scale.
#'
#' @details
#' This function implements the high-dimensional two-step estimator for linear
#' expected shortfall regression proposed by Zhang et al. (2025). The input
#' matrix \code{x} should not contain an intercept column. Intercepts for the
#' quantile and expected shortfall regressions are estimated separately.
#'
#' In the first step, a penalized quantile regression is fitted at level
#' \code{alpha}. If \code{tuning.method = "CV"}, the tuning parameter is
#' selected by cross-validation using
#' \code{\link[conquer]{conquer.cv.reg}}. If
#' \code{tuning.method = "BIC"}, the tuning parameter is selected using a
#' BIC-type criterion based on \code{\link[conquer]{conquer.reg}}.
#'
#' Let \eqn{\hat{Q}_i} denote the fitted conditional quantile for the \eqn{i}-th
#' observation. The pseudo-response
#' \deqn{
#' \hat{Z}_i = \hat {Q}_i + \alpha^{-1} \mathbf{1}(Y_i \leq \hat{Q}_i)(Y_i - \hat{Q}_i)
#' }
#' is then regressed on the same design matrix by penalized least squares using
#' \code{\link[glmnet]{glmnet}}. The resulting coefficients estimate the
#' high-dimensional linear expected shortfall function.
#'
#' If \code{upperTail = TRUE}, the upper-tail model is fitted using the same
#' transformation principle as in \code{\link{lesreg}}. Statistical inference
#' for an individual ES regression coefficient can be obtained using
#' \code{\link{lesreg.highdim.inference}}.
#'
#' @return
#' \code{lesreg.highdim} returns an object of class
#' \code{"lesreg.highdim"}.
#'
#' An object of class \code{"lesreg.highdim"} is a list containing at least the
#' following components:
#' \describe{
#'   \item{\code{theta0_hat}}{the intercept from the second-step expected
#'   shortfall regression fit.}
#'   \item{\code{theta_hat}}{the slope coefficients from the second-step
#'   expected shortfall regression fit.}
#'   \item{\code{beta0_hat}}{the intercept from the first-step quantile
#'   regression fit.}
#'   \item{\code{beta_hat}}{the slope coefficients from the first-step quantile
#'   regression fit.}
#'   \item{\code{resid}}{the residuals from the second-step penalized
#'   least-squares fit, defined as the pseudo-response minus the fitted expected
#'   shortfall values.}
#'   \item{\code{lambda_beta}}{the selected tuning parameter from the
#'   first-step quantile regression fit.}
#'   \item{\code{lambda_theta}}{the selected tuning parameter from the
#'   second-step expected shortfall regression fit.}
#'   \item{\code{tau}}{the tail probability used for estimation.}
#'   \item{\code{upperTail}}{logical flag indicating whether the upper tail was
#'   fitted.}
#'   \item{\code{tail}}{a character string indicating whether the lower or
#'   upper tail was fitted.}
#'   \item{\code{standardize}}{logical flag indicating whether standardization
#'   was requested for the penalized fitting procedures.}
#'   \item{\code{tuning.method}}{the tuning parameter selection method used.}
#'   \item{\code{call}}{the matched call.}
#'   \item{\code{y}}{the original response vector.}
#'   \item{\code{y.fit}}{the response vector used internally for fitting. This
#'   is \code{y} for lower-tail fitting and \code{-y} for upper-tail fitting.}
#'   \item{\code{theta_hat_fit}}{the second-step slope coefficients on the
#'   internal fitting scale.}
#'   \item{\code{theta0_hat_fit}}{the second-step intercept on the internal
#'   fitting scale.}
#'   \item{\code{beta_hat_fit}}{the first-step slope coefficients on the
#'   internal fitting scale.}
#'   \item{\code{beta0_hat_fit}}{the first-step intercept on the internal
#'   fitting scale.}
#'   \item{\code{resid_fit}}{the second-step residuals on the internal fitting
#'   scale.}
#' }
#'
#' @references
#' Zhang, S., He, X., Tan, K. M., and Zhou, W.-X. (2025).
#' High-Dimensional Expected Shortfall Regression.
#' \emph{Journal of the American Statistical Association}, 120(551), 1799--1810.
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
#' ## Fit the model with l1 penalties selected by cross-validation in two steps
#' fit <- lesreg.highdim(x, y, alpha = 0.05)
#' summary(fit)
#'
#' ## Fit the model with l1 penalties selected by BIC in two steps
#' fit.bic <- lesreg.highdim(x, y, alpha = 0.05, tuning.method = "BIC")
#' summary(fit.bic)
#'
#' ## Fit the model for upper-tail ES
#' fit.upper <- lesreg.highdim(x, y, alpha = 0.05, upperTail = TRUE)
#' summary(fit.upper)
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
