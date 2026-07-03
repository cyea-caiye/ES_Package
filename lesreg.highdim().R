#' @export
lesreg.highdim <- function(x, y, alpha = 0.05, standardize = TRUE, tuning.method = c("CV", "BIC")){
  call <- match.call()
  tuning.method <- match.arg(tuning.method)
  if(is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  p <- ncol(x)
  if(tuning.method == "BIC"){
    quant_model_bic <- .lesreg_hd_quantile_bic(
      x,
      y,
      alpha,
      standardize = standardize
    )
    lambda_beta <- quant_model_bic$lambda_select
    model_select <- quant_model_bic$model_select
    beta_hat <- model_select$coeff[-1]
    beta0 <- model_select$coeff[1]
    Z_2step <- apply(
      cbind(y, x),
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
      y,
      tau = alpha,
      h = max(0.05, sqrt(alpha * (1 - alpha)) * (log(p) / n)^(1 / 4))
    )
    lambda_beta <- quan_model$lambda
    beta_hat <- quan_model$coeff.min[-1]
    beta0 <- quan_model$coeff.min[1]
    Z_2step <- beta0 + x %*% beta_hat +
      (y - beta0 - x %*% beta_hat) *
      ifelse(y < beta0 + x %*% beta_hat, 1, 0) / alpha
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
  ans <- list(
    theta0_hat = theta0_hat,
    theta_hat = theta_hat,
    beta0_hat = beta0,
    beta_hat = beta_hat,
    resid = resid,
    lambda_beta = lambda_beta,
    lambda_theta = lambda_theta,
    alpha = alpha,
    standardize = standardize,
    tuning.method = tuning.method,
    call = call
  )
  class(ans) <- "lesreg.highdim"
  ans
}