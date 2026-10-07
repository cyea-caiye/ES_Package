#' Summarizing High-Dimensional Linear Expected Shortfall Regression Fits
#'
#' @description
#' Summarizes a fitted high-dimensional linear expected shortfall (ES)
#' regression model returned by \code{\link{lesreg.highdim}}. The summary
#' reports the fitted quantile and ES regression coefficients, selected tuning
#' parameters, and the numbers of nonzero slope coefficients.
#'
#' @param object an object of class \code{"lesreg.highdim"}, usually a result
#'   of a call to \code{\link{lesreg.highdim}}.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The fitted object contains penalized coefficient estimates from both steps
#' of the high-dimensional ES regression procedure. This method organizes the
#' first-step quantile regression coefficients and second-step ES regression
#' coefficients, including their respective intercepts, and reports the number
#' of selected nonzero slope coefficients in each step.
#'
#' See \code{\link{lesreg.highdim}} for details on the underlying two-step
#' estimation procedure and tuning parameter selection.
#'
#' @return
#' \code{summary.lesreg.highdim} returns an object of class
#' \code{"summary.lesreg.highdim"} containing:
#' \describe{
#'   \item{\code{call}}{the matched call from the fitted model.}
#'   \item{\code{tau}}{the tail probability used for estimation.}
#'   \item{\code{upperTail}}{logical flag indicating whether the upper tail was
#'   fitted.}
#'   \item{\code{tuning.method}}{the tuning parameter selection method used.}
#'   \item{\code{standardize}}{logical flag indicating whether standardization
#'   was requested for the penalized fitting procedures.}
#'   \item{\code{lambda_beta}}{the selected tuning parameter from the
#'   first-step quantile regression fit.}
#'   \item{\code{lambda_theta}}{the selected tuning parameter from the
#'   second-step expected shortfall regression fit.}
#'   \item{\code{quantile.coefficients}}{a named vector containing the
#'   first-step quantile regression intercept and slope coefficients.}
#'   \item{\code{es.coefficients}}{a named vector containing the second-step
#'   expected shortfall regression intercept and slope coefficients.}
#'   \item{\code{nzero_beta}}{the number of nonzero slope coefficients in the
#'   first-step quantile regression fit.}
#'   \item{\code{nzero_theta}}{the number of nonzero slope coefficients in the
#'   second-step expected shortfall regression fit.}
#'   \item{\code{p}}{the number of slope coefficients in the fitted model.}
#' }
#'
#' @seealso
#' \code{\link{lesreg.highdim}},
#' \code{\link{lesreg.highdim.inference}}
#'
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
    upperTail = object$upperTail,
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
