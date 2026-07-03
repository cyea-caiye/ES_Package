#' @export
lesreg.highdim.inf <- function(x, y, alpha, col = 1, res_est = NULL, conf.level = 0.95, standardize = TRUE,
                               variance.method = c("RCV", "refit", "plug-in")){
  variance.method <- match.arg(variance.method)
  if(is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  p <- ncol(x)
  if(col < 1 || col > p || col != as.integer(col))
    stop("'col' must be an integer between 1 and ncol(x)")
  if(!is.numeric(conf.level) || length(conf.level) != 1L ||
     is.na(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("'conf.level' must be a single number between 0 and 1")
  if(!is.logical(standardize) || length(standardize) != 1L ||
     is.na(standardize))
    stop("'standardize' must be TRUE or FALSE")
  if(is.null(res_est)){
    res_est <- lesreg.highdim(x, y, alpha = alpha, standardize = standardize)
  }
  theta0_hat <- res_est$theta0_hat
  theta_hat <- res_est$theta_hat
  resid1 <- res_est$resid
  beta0 <- res_est$beta0_hat
  beta_hat <- res_est$beta_hat
  d <- x[, col]
  x_tilde <- x[, -col, drop = FALSE]
  ### CV ###
  cv_model <- glmnet::cv.glmnet(x_tilde, d, standardize = standardize)
  lambda_cv <- cv_model$lambda.1se
  gamma_model <- glmnet::glmnet(x_tilde, d, lambda = lambda_cv, standardize = standardize)
  gamma_hat <- gamma_model$beta
  nonzero_index <- which(gamma_hat != 0)
  resid2 <- d - stats::predict(gamma_model, newx = x_tilde)
  w_hat <- resid2
  theta_debias <- theta_hat[col] + t(resid1) %*% resid2 / x[, col] %*% resid2
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
    Z_2step <- apply(cbind(y, x), MARGIN = 1, FUN = .lesreg_hd_Z_beta, beta = beta_hat, alpha = alpha, beta_0 = beta0)
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
    y_1 <- y[sample_index]
    y_2 <- y[-sample_index]
    var <- .lesreg_hd_variance_estimation(d_x_1, d_x_2, y_1, y_2, alpha, standardize = standardize)
    wid <- sqrt(var) * stats::qnorm(CI_level) / sqrt(n) / alpha
  }
  return(list(
    theta_debias = theta_debias,
    Conf.int = c(theta_debias - wid, theta_debias + wid)
  ))
}
