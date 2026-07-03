.lesreg_hd_Z_beta <- function(data, beta, alpha, beta_0){
  yi <- data[1]
  xi <- data[-1]
  res <- yi - xi %*% beta - beta_0
  if(res <= 0){
    return(res / alpha + xi %*% beta + beta_0)
  } else {
    return(xi %*% beta + beta_0)
  }
}

.lesreg_hd_self_tuning <- function(x, tau, nlambda, standardize = TRUE){
  n <- nrow(x)
  if(standardize){
    x <- scale(x, center = TRUE, scale = TRUE)
  }
  lambda_sim <- numeric(nlambda)
  for(b in 1:nlambda){
    u <- stats::runif(nrow(x), 0, 1) <= tau
    lambda_sim[b] <- max(abs(crossprod(x, tau - u)))
  }
  return(lambda_sim / n)
}

.lesreg_hd_quantile_bic <- function(x, y, tau, lambda_seq = numeric(), nlambda = 100, h = NULL, Cn = NULL, standardize = TRUE){
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(h)){
    h <- max(0.05, sqrt(tau * (1 - tau)) * (log(p) / n)^(1 / 4))
  }
  if(is.null(Cn)){
    Cn <- max(2, log(log(n)))
  }
  if(length(lambda_seq) == 0){
    lam_max <- max(.lesreg_hd_self_tuning(x, tau, nlambda, standardize = standardize))
    lambda_seq <- seq(0.25 * lam_max, lam_max, length.out = nlambda)
  } else {
    nlambda <- length(lambda_seq)
  }
  check_sum <- function(x, tau){
    sum(ifelse(x >= 0, tau * x, (tau - 1) * x))
  }
  model_all <- conquer::conquer.reg(x, y, lambda = lambda_seq, tau = tau, h = h)
  beta_hat <- model_all$coeff
  BIC <- log(apply(matrix(rep(y, nlambda), ncol = nlambda) - matrix(rep(beta_hat[1, ], each = n), ncol = nlambda) - x %*% beta_hat[-1, ], 2, check_sum, tau = tau)) + apply(beta_hat, 2, function(x){sum(x != 0)}) * log(p) * Cn / n
  bic_select <- which.min(BIC)
  lambda_select <- lambda_seq[bic_select]
  model_select <- conquer::conquer.reg(x, y, lambda = lambda_select, tau = tau, h = h)
  return(list(lambda_select = lambda_select, model_select = model_select))
}

.lesreg_hd_LS_bic <- function(x, z, lambda_seq = numeric(), nlambda = 100, Cn = NULL, epsilon = .0001, standardize = TRUE){
  n <- nrow(x)
  p <- ncol(x)
  if(is.null(Cn)){
    Cn <- max(2, log(log(n)))
  }
  model_all <- glmnet::glmnet(x, z, nlambda = nlambda, standardize = standardize)
  theta_hat <- model_all$beta
  theta0 <- model_all$a0
  lambda_seq <- model_all$lambda
  nlambda <- length(lambda_seq)
  BIC <- log(apply((matrix(rep(z, nlambda), ncol = nlambda) - matrix(rep(theta0, each = n), ncol = nlambda) - x %*% theta_hat)^2, 2, sum) / 2 / n) + apply(theta_hat, 2, function(x){sum(x != 0)}) * log(p) * Cn / n
  bic_select <- which.min(BIC)
  lambda_select <- lambda_seq[bic_select]
  model_select <- glmnet::glmnet(x, z, lambda = lambda_select, standardize = standardize)
  return(list(lambda_select = lambda_select, model_select = model_select))
}

.lesreg_hd_variance_estimation <- function(d_x_1, d_x_2, y_1, y_2, alpha, standardize = TRUE){
  n <- nrow(d_x_1) + nrow(d_x_2)
  p <- ncol(d_x_1)
  # Remove the columns that are constant
  constant_columns_1 = data.frame(d_x_1) %>%
    dplyr::select_at(setdiff(names(.), names(janitor::remove_constant(.)))) %>%
    unique()
  constant_columns_2 = data.frame(d_x_2) %>%
    dplyr::select_at(setdiff(names(.), names(janitor::remove_constant(.)))) %>%
    unique()
  constant_col = unique(c(colnames(constant_columns_1), colnames(constant_columns_2)))
  d_x_1 = as.matrix(dplyr::select(data.frame(d_x_1), -constant_col))
  d_x_2 = as.matrix(dplyr::select(data.frame(d_x_2), -constant_col))
  #######find variables on first half######## 
  #######three step procedure#########
  # step 1
  quan_model = conquer::conquer.cv.reg(d_x_1, y_1,
                                       #lambda=0.03,
                                       tau = alpha,
                                       h = max(0.05, sqrt(alpha * (1 - alpha)) * (log(p) / n)^(1 / 4)))
  beta_hat = quan_model$coeff.min[-1]
  beta0 = quan_model$coeff.min[1]
  sq = sum(beta_hat != 0)
  # step 2
  Z_2step = apply(cbind(y_1, d_x_1), MARGIN = 1, FUN = .lesreg_hd_Z_beta, beta = beta_hat, alpha = alpha, beta_0 = beta0)
  ######### CV #########
  ##perform k-fold cross-validation to find optimal lambda value
  cv_model <- glmnet::cv.glmnet(d_x_1, Z_2step, standardize = standardize)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet::glmnet(d_x_1, Z_2step, lambda = lambda_theta, standardize = standardize)
  theta_hat = theta_model$beta
  se = length(theta_hat@i)
  # step 3
  #d_x_1_stand = sweep(d_x_1, 2, colMeans(d_x_1))
  d_x_1_stand = d_x_1
  cv_model <- glmnet::cv.glmnet(d_x_1_stand[, -1], d_x_1_stand[, 1], standardize = standardize)
  lambda_cv <- cv_model$lambda.1se
  dx1_model = glmnet::glmnet(d_x_1_stand[, -1], d_x_1_stand[, 1], lambda = lambda_cv, standardize = standardize)
  gamma_hat = dx1_model$beta
  nonzero_index = which(gamma_hat != 0)
  #print(nonzero_index)
  #nonzero_index = gamma_hat@i+2
  sm = length(nonzero_index)
  ##########fit on second half###################
  # step 1
  if(sq == 0){
    epsilon_hat = y_2 - beta0
    Z_2step = beta0 + ifelse(epsilon_hat < 0, epsilon_hat, 0) / alpha
  } else{
    dx2_short = matrix(d_x_2[, beta_hat != 0], ncol = sq)
    quan_model <- conquer::conquer(dx2_short, y_2, tau = alpha)
    beta_hat = quan_model$coeff[-1]
    beta0 = quan_model$coeff[1]
    #epsilon_hat = y_2 - dx2_short%*%beta_hat - beta0
    #epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
    Z_2step = apply(cbind(y_2, dx2_short), MARGIN = 1, FUN = .lesreg_hd_Z_beta, beta = beta_hat, alpha = alpha, beta_0 = beta0)
  }
  # step 2
  if(se == 0){
    res3 = Z_2step - mean(Z_2step)
  } else{
    dx2_short_theta = matrix(d_x_2[, theta_hat@i + 1], ncol = se)
    theta_model = stats::lm(Z_2step ~ dx2_short_theta)
    theta_hat = theta_model$coefficients[-1]
    res3 = theta_model$residuals
  }
  # step 3
  #d_x_2_stand = sweep(d_x_2, 2, colMeans(d_x_2))
  d_x_2_stand = d_x_2
  if(sm == 0){
    w_hat = d_x_2_stand[, 1] - mean(d_x_2_stand[, 1])
  } else{
    lm_model <- stats::lm(d_x_2_stand[, 1] ~ d_x_2_stand[, nonzero_index + 1])
    w_hat = lm_model$residuals
  }
  #results
  s2 = nrow(d_x_2)
  sigma_w_hat_square_1 = sum(w_hat^2) / (s2 - sm - 1)
  sigma_s_hat_1 = alpha^2 * t(w_hat^2) %*% ((res3)^2) / (s2 - sm - sq - se - 3)
  
  #######################The other way around#####################
  #######find variables on first half######## 
  #######three step procedure#########
  
  # step 1
  quan_model = conquer::conquer.cv.reg(d_x_2, y_2,
                                       #lambda=0.03,
                                       tau = alpha,
                                       h = max(0.05, sqrt(alpha * (1 - alpha)) * (log(p) / n)^(1 / 4)))
  beta_hat = quan_model$coeff.min[-1]
  beta0 = quan_model$coeff.min[1]
  sq = sum(beta_hat != 0)
  # step 2
  Z_2step = apply(cbind(y_2, d_x_2), MARGIN = 1, FUN = .lesreg_hd_Z_beta, beta = beta_hat, alpha = alpha, beta_0 = beta0)
  ######### CV #########
  ##perform k-fold cross-validation to find optimal lambda value
  cv_model <- glmnet::cv.glmnet(d_x_2, Z_2step, standardize = standardize)
  ##find optimal lambda value that minimizes test MSE
  lambda_theta <- cv_model$lambda.min
  theta_model = glmnet::glmnet(d_x_2, Z_2step, lambda = lambda_theta, standardize = standardize)
  theta_hat = theta_model$beta
  se = length(theta_hat@i)
  # step 3
  cv_model <- glmnet::cv.glmnet(d_x_2_stand[, -1], d_x_2_stand[, 1], standardize = standardize)
  lambda_cv <- cv_model$lambda.1se
  dx2_model = glmnet::glmnet(d_x_2_stand[, -1], d_x_2_stand[, 1], lambda = lambda_cv, standardize = standardize)
  gamma_hat = dx2_model$beta
  nonzero_index = which(gamma_hat != 0)
  #print(nonzero_index)
  #nonzero_index = gamma_hat@i+2
  sm = length(nonzero_index)
  ##########fit on second half###################
  # step 1
  if(sq == 0){
    epsilon_hat = y_1 - beta0
    Z_2step = beta0 + ifelse(epsilon_hat < 0, epsilon_hat, 0) / alpha
  } else{
    dx1_short = matrix(d_x_1[, beta_hat != 0], ncol = sq)
    quan_model <- conquer::conquer(dx1_short, y_1, tau = alpha)
    beta_hat = quan_model$coeff[-1]
    beta0 = quan_model$coeff[1]
    #epsilon_hat = y_1 - dx1_short%*%beta_hat - beta0
    #epsilon_minus = ifelse(epsilon_hat<0,epsilon_hat,0)
    Z_2step = apply(cbind(y_1, dx1_short), MARGIN = 1, FUN = .lesreg_hd_Z_beta, beta = beta_hat, alpha = alpha, beta_0 = beta0)
  }
  # step 2
  if(se == 0){
    res3 = Z_2step - mean(Z_2step)
  } else{
    dx1_short_theta = matrix(d_x_1[, theta_hat@i + 1], ncol = se)
    theta_model = stats::lm(Z_2step ~ dx1_short_theta)
    theta_hat = theta_model$coefficients[-1]
    res3 = theta_model$residuals
  }
  # step 3
  if(sm == 0){
    w_hat = d_x_1_stand[, 1] - mean(d_x_1_stand[, 1])
  } else{
    lm_model <- stats::lm(d_x_1_stand[, 1] ~ d_x_1_stand[, nonzero_index + 1])
    w_hat = lm_model$residuals
  }
  #results
  s1 = nrow(d_x_1)
  sigma_w_hat_square_2 = sum(w_hat^2) / (s1 - sm - 1)
  sigma_s_hat_2 = alpha^2 * t(w_hat^2) %*% ((res3)^2) / (s1 - sm - sq - se - 3)
  sigma_w_hat_square = (sigma_w_hat_square_1 + sigma_w_hat_square_2) / 2
  sigma_s_hat = (sigma_s_hat_1 + sigma_s_hat_2) / 2
  #CI_level = 1-(1-conf.level)/2
  var = sigma_s_hat / sigma_w_hat_square^2
  return(var)
}