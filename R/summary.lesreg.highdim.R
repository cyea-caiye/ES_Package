#' @rdname lesreg.highdim
#' @param object an object of class \code{"lesreg.highdim"}.
#' @export
summary.lesreg.highdim <- function(object, ...){
  get_coef <- function(intercept, beta){
    if(inherits(beta, "Matrix")){
      beta_names <- rownames(beta)
      beta <- as.vector(beta)
      names(beta) <- beta_names
    } else{
      beta <- as.vector(beta)
    }
    if(is.null(names(beta))){
      names(beta) <- paste0("V", seq_along(beta))
    }
    c("(Intercept)" = as.numeric(intercept), beta)
  }
  beta_coef <- get_coef(object$beta0_hat, object$beta_hat)
  theta_coef <- get_coef(object$theta0_hat, object$theta_hat)
  ans <- list(
    call = object$call,
    tau = object$tau,
    uppertail = object$uppertail,
    tuning.method = object$tuning.method,
    standardize = object$standardize,
    lambda_beta = object$lambda_beta,
    lambda_theta = object$lambda_theta,
    quantile.coefficients = beta_coef,
    es.coefficients = theta_coef,
    nzero_beta = sum(beta_coef[-1] != 0),
    nzero_theta = sum(theta_coef[-1] != 0),
    p = length(beta_coef) - 1
  )
  class(ans) <- "summary.lesreg.highdim"
  ans
}
