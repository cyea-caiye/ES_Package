#' Inference for High-Dimensional Linear Expected Shortfall Regression Models
#'
#' @description
#' \code{lesreg.highdim.inference} performs inference on one coefficient of a
#' high-dimensional linear expected shortfall regression model using a debiased
#' estimator.
#'
#' @param x numeric design matrix of dimension \eqn{n * p}, without an
#'   intercept column.
#' @param y vector of observations of length \eqn{n}, or a one-column matrix.
#' @param alpha the tail probability used for expected shortfall regression. It
#'   must be a single number strictly between 0 and 1. If
#'   \code{upperTail = FALSE}, \code{alpha} is the lower-tail probability. If
#'   \code{upperTail = TRUE}, \code{alpha} is the upper-tail probability.
#' @param col integer specifying the index of the coefficient of interest.
#' @param res_est optional fitted object returned by \code{\link{lesreg.highdim}}.
#'   If \code{NULL}, \code{lesreg.highdim} is called internally.
#' @param conf.level confidence level for the reported confidence interval.
#' @param standardize logical. If \code{TRUE}, covariates are standardized in the
#'   penalized fitting routines.
#' @param variance.method character string specifying the method used to
#'   estimate the asymptotic variance. Available choices are \code{"RCV"},
#'   \code{"refit"} and \code{"plug-in"}.
#' @param upperTail logical. If \code{FALSE}, inference is performed for the
#'   lower-tail expected shortfall coefficient. If \code{TRUE}, inference is
#'   performed for the upper-tail coefficient by applying the lower-tail
#'   procedure to the negative response and flipping the sign of the debiased
#'   estimator and confidence interval.
#'
#' @details
#' The fitted high-dimensional expected shortfall regression coefficients are
#' penalized estimators and are not directly used for ordinary low-dimensional
#' inference. For a specified coefficient indexed by \code{col}, this function
#' computes a debiased estimator by projecting the covariate of interest on the
#' remaining covariates and using the resulting residuals to form a
#' decorrelated score.
#'
#' If \code{res_est} is supplied, the first-step quantile regression and
#' second-step expected shortfall regression estimates are taken from that
#' object. Otherwise, \code{\link{lesreg.highdim}} is called internally.
#'
#' Three variance estimation methods are available. The method \code{"plug-in"}
#' uses a direct plug-in variance estimate. The method \code{"refit"} refits
#' selected low-dimensional models on the selected supports. The method
#' \code{"RCV"} uses a refitted cross-validation variance estimator.
#'
#' @return
#' \code{lesreg.highdim.inference} returns an object of class
#' \code{"lesreg.highdim.inference"}.
#'
#' An object of class \code{"lesreg.highdim.inference"} is a list containing
#' \describe{
#'   \item{theta_debias}{the debiased estimator for the coefficient of
#'   interest.}
#'   \item{Conf.int}{confidence interval for the debiased estimator.}
#'   \item{col}{the index of the coefficient of interest.}
#'   \item{tau}{the tail probability used for estimation.}
#'   \item{conf.level}{the confidence level used for the reported confidence
#'   interval.}
#'   \item{variance.method}{the variance estimation method used.}
#'   \item{call}{the matched call.}
#'   \item{upperTail}{logical flag indicating whether inference was performed
#'   for the upper tail.}
#'   \item{tail}{character string indicating whether inference was performed for
#'   the \code{"lower"} or \code{"upper"} tail.}
#' }
#'
#' @seealso
#' \code{\link{lesreg.highdim}}, \code{\link{lesreg}},
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
#'
#' inf <- lesreg.highdim.inference(
#'   x, y,
#'   alpha = 0.05,
#'   col = 1,
#'   res_est = fit,
#'   variance.method = "plug-in"
#' )
#' inf
#'
#' fit.upper <- lesreg.highdim(x, y, alpha = 0.05, upperTail = TRUE)
#' inf.upper <- lesreg.highdim.inference(
#'   x, y,
#'   alpha = 0.05,
#'   col = 1,
#'   res_est = fit.upper,
#'   variance.method = "plug-in"
#' )
#' inf.upper
#'
#' @export
lesreg.highdim.inference <- function(x, y, alpha, col = 1, res_est = NULL, conf.level = 0.95, standardize = TRUE,
                               variance.method = c("RCV", "refit", "plug-in"), upperTail = FALSE){
  variance.method <- match.arg(variance.method)
  if(is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  p <- ncol(x)
  if(col < 1 || col > p || col != as.integer(col))
    stop("'col' must be an integer between 1 and ncol(x)")
  if(!is.numeric(conf.level) || length(conf.level) != 1L ||
     is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number between 0 and 1")
  if(!is.logical(standardize) || length(standardize) != 1L || is.na(standardize))
    stop("'standardize' must be TRUE or FALSE")
  if(!is.logical(upperTail) || length(upperTail) != 1L || is.na(upperTail)){
    stop("'upperTail' must be TRUE or FALSE")
  }
  if(is.null(res_est)){
    res_est <- lesreg.highdim(x, y, alpha = alpha, standardize = standardize, upperTail = upperTail)
  } else{
    upperTail <- res_est$upperTail
  }
  y_fit <- res_est$y.fit
  theta0_hat <- res_est$theta0_hat_fit
  theta_hat <- res_est$theta_hat_fit
  resid1 <- res_est$resid_fit
  beta0 <- res_est$beta0_hat_fit
  beta_hat <- res_est$beta_hat_fit
  d <- x[, col]
  x_tilde <- x[, -col, drop = FALSE]
  ### CV ###
  cv_model <- glmnet::cv.glmnet(x_tilde, d, standardize = standardize)
  lambda_cv <- cv_model$lambda.1se
  gamma_model <- glmnet::glmnet(x_tilde, d, lambda = lambda_cv, standardize = standardize)
  gamma_hat <- gamma_model$beta
  nonzero_index = which(as.vector(gamma_hat) != 0)
  resid2 <- d - stats::predict(gamma_model, newx = x_tilde)
  w_hat <- resid2
  theta_debias_fit <- theta_hat[col] + t(resid1) %*% resid2 / x[, col] %*% resid2
  sq <- sum(beta_hat != 0)
  se <- length(theta_hat@i)
  sm <- length(nonzero_index)
  CI_level <- 1 - (1 - conf.level) / 2
  ############ Variance estimation using simple plug-in ################
  if(variance.method == "plug-in"){
    sigma_w_hat_square <- sum(w_hat^2) / (n - sm - 1)
    sigma_s_hat <- alpha^2 * t(w_hat^2) %*% ((resid1)^2) / (n - sm - sq - se - 3)
    var_plugin <- sigma_s_hat / sigma_w_hat_square^2
    wid <- sqrt(var_plugin) * stats::qnorm(CI_level) / sqrt(n) / alpha
  } else if(variance.method == "refit"){
    ################ Variance Estimation using vanilla refit ################
    Z_2step <- apply(cbind(y_fit, x), MARGIN = 1, FUN = .lesreg_hd_Z_beta, beta = beta_hat, alpha = alpha, beta_0 = beta0)
    if(se == 0){
      res3 <- Z_2step - mean(Z_2step)
    } else{
      x_refit <- x[, theta_hat@i + 1]
      theta_model_refit <- stats::lm(Z_2step ~ x_refit)
      theta_hat_refit <- theta_model_refit$coefficients[-1]
      res3 <- theta_model_refit$residuals
    }
    if(sm == 0){
      w_hat <- d - mean(d)
    } else{
      lm_model <- stats::lm(d ~ x_tilde[, nonzero_index])
      w_hat <- lm_model$residuals
    }
    sigma_w_hat_square <- sum(w_hat^2) / (n - sm - 1)
    sigma_s_hat <- alpha^2 * t(w_hat^2) %*% ((res3)^2) / (n - sm - sq - se - 3)
    var_vanilla <- sigma_s_hat / sigma_w_hat_square^2
    wid <- sqrt(var_vanilla) * stats::qnorm(CI_level) / sqrt(n) / alpha
  } else{
    ################ Variance Estimation using RCV ################
    sample_index <- sample(seq_len(nrow(x)), size = ceiling(n / 2))
    d_x_1 <- x[sample_index, ]
    d_x_2 <- x[-sample_index, ]
    y_1 <- y_fit[sample_index]
    y_2 <- y_fit[-sample_index]
    var <- .lesreg_hd_variance_estimation(d_x_1, d_x_2, y_1, y_2, alpha, col, standardize = standardize)
    wid <- sqrt(var) * stats::qnorm(CI_level) / sqrt(n) / alpha
  }
  theta_debias_fit <- as.numeric(theta_debias_fit)
  wid <- as.numeric(wid)
  if(upperTail){
    theta_debias <- -theta_debias_fit
  } else {
    theta_debias <- theta_debias_fit
  }
  ans <- list(
    theta_debias = theta_debias,
    Conf.int = c(theta_debias - wid, theta_debias + wid),
    col = col,
    tau = alpha,
    upperTail = upperTail,
    tail = if(upperTail) "upper" else "lower",
    conf.level = conf.level,
    variance.method = variance.method,
    call = match.call()
  )
  class(ans) <- "lesreg.highdim.inference"
  ans
}
