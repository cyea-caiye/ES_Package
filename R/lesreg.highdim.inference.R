#' Inference for High-Dimensional Linear Expected Shortfall Regression Models
#'
#' @description
#' Performs statistical inference for an individual coefficient in a
#' high-dimensional linear expected shortfall (ES) regression model using a
#' debiased estimator. The procedure is based on the high-dimensional ES
#' regression framework of Zhang et al. (2025).
#'
#' @param x numeric design matrix of dimension \eqn{n \times p}, without an
#'   intercept column.
#' @param y numeric response vector of length \eqn{n}, or a one-column numeric
#'   matrix.
#' @param alpha the quantile level associated with the ES to be estimated. This
#'   must be a single number strictly between 0 and 1. If
#'   \code{upperTail = FALSE}, the function performs inference for the
#'   lower-tail ES associated with quantile level \code{alpha}. If
#'   \code{upperTail = TRUE}, the function performs inference for the
#'   upper-tail ES associated with quantile level \code{1 - alpha}. This
#'   corresponds to \code{tau} in \code{\link{lesreg}}.
#' @param col integer specifying the column of \code{x}, and hence the slope
#'   coefficient, for which inference is performed. This must be an integer
#'   between 1 and \code{ncol(x)}.
#' @param res_est optional fitted object returned by
#'   \code{\link{lesreg.highdim}}. If \code{NULL},
#'   \code{lesreg.highdim} is called internally.
#' @param conf.level confidence level for the reported confidence interval.
#'   This must be a single number strictly between 0 and 1.
#' @param standardize logical. If \code{TRUE}, standardization is used in the
#'   penalized fitting procedures where supported.
#' @param variance.method character string specifying the method used to
#'   estimate the asymptotic variance. Available choices are \code{"RCV"},
#'   \code{"refit"}, and \code{"plug-in"}.
#' @param upperTail logical. If \code{FALSE}, inference is performed for the
#'   lower-tail ES coefficient. If \code{TRUE}, inference is performed for the
#'   upper-tail ES coefficient by applying the lower-tail procedure to the
#'   negative response and transforming the resulting debiased estimator back
#'   to the original response scale.
#'
#' @details
#' The penalized coefficient estimates returned by
#' \code{\link{lesreg.highdim}} are regularized and are therefore not directly
#' used for conventional low-dimensional inference. For the slope coefficient
#' indexed by \code{col}, \code{lesreg.highdim.inference} constructs a
#' debiased estimator using a decorrelated-score approach. Specifically, the
#' covariate of interest is projected onto the remaining covariates using a
#' penalized regression, and the resulting projection residuals are used to
#' correct the initial penalized ES coefficient estimator.
#'
#' If \code{res_est} is supplied, the first-step quantile regression and
#' second-step expected shortfall regression estimates are taken from that
#' fitted object. Otherwise, \code{\link{lesreg.highdim}} is called internally
#' using \code{x}, \code{y}, \code{alpha}, \code{standardize}, and
#' \code{upperTail}.
#'
#' Three variance estimation methods are available. The
#' \code{"plug-in"} method uses a direct plug-in estimator of the asymptotic
#' variance. The \code{"refit"} method refits low-dimensional regression models
#' on the supports selected by the penalized estimators before estimating the
#' variance. The \code{"RCV"} method uses a refitted cross-validation variance
#' estimator based on a random sample split.
#'
#' If \code{upperTail = TRUE}, the inference procedure is carried out on the
#' internal lower-tail representation based on the negative response, and the
#' debiased estimator is transformed back to the original response scale. See
#' \code{\link{lesreg.highdim}} and \code{\link{lesreg}} for further details on
#' the upper-tail transformation.
#'
#' @return
#' \code{lesreg.highdim.inference} returns an object of class
#' \code{"lesreg.highdim.inference"}.
#'
#' An object of class \code{"lesreg.highdim.inference"} is a list containing
#' the following components:
#' \describe{
#'   \item{\code{theta_debias}}{the debiased estimator for the slope coefficient
#'   indexed by \code{col}.}
#'   \item{\code{Conf.int}}{the confidence interval for the debiased estimator.}
#'   \item{\code{col}}{the index of the coefficient of interest.}
#'   \item{\code{tau}}{the tail probability used for estimation.}
#'   \item{\code{conf.level}}{the confidence level used for the reported
#'   confidence interval.}
#'   \item{\code{variance.method}}{the variance estimation method used.}
#'   \item{\code{call}}{the matched call.}
#'   \item{\code{upperTail}}{logical flag indicating whether inference was
#'   performed for the upper tail.}
#'   \item{\code{tail}}{a character string indicating whether inference was
#'   performed for the lower or upper tail.}
#' }
#'
#' @references
#' Zhang, S., He, X., Tan, K. M., and Zhou, W.-X. (2025).
#' High-Dimensional Expected Shortfall Regression.
#' \emph{Journal of the American Statistical Association}, 120(551), 1799--1810.
#'
#' @seealso
#' \code{\link{lesreg.highdim}}, \code{\link{lesreg}},
#' \code{\link[glmnet]{glmnet}},
#' \code{\link[glmnet]{cv.glmnet}}
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' p <- 100
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - x[, 2] + rnorm(n)
#' fit <- lesreg.highdim(x, y, alpha = 0.05)
#'
#' ## Conduct statistical inference on the 1st, 2nd, 3rd, 5th and 7th coefficients
#' coefs.true <- c(1,-1, 0,0,0)
#' coef.debiase <- c()
#' cis <- NULL
#'
#' for(col in c(1,2,3,5,7)) {
#'   inf <- lesreg.highdim.inference(
#'     x, y,
#'     alpha = 0.05,
#'     col = col,
#'     res_est = fit,
#'     variance.method = "RCV"
#'   )
#'   coef.debiase <- c(coef.debiase, inf$theta_debias)
#'   cis <- cbind(cis, inf$Conf.int)
#' }
#'
#' res.table <- rbind(coefs.true, coef.debiase, cis)
#' rownames(res.table) <- c('coef.true', 'coef.debiase', 'ci.low', 'ci.up')
#' print(res.table)
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
  theta_debias_fit <- theta_hat[col] + t(resid1) %*% resid2 / (x[, col] %*% resid2)
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
