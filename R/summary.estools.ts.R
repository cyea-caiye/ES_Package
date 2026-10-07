#' Summarizing Autoregressive Quantile-Expected Shortfall Fits
#'
#' @description
#' \code{summary.estools.ts} computes inference summaries for the final
#' weighted expected shortfall estimator returned by
#' \code{\link{estools.ts}}.
#'
#' @param object an object of class \code{"estools.ts"}.
#' @param level confidence level for the reported confidence intervals.
#'   It must be a single number strictly between 0 and 1.
#' @param covariance logical. If \code{TRUE}, the estimated coefficient
#'   covariance matrix together with the sandwich components
#'   \code{Sigma_hat} and \code{Omega_hat} is included in the returned object.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' Inference is reported for the final weighted expected shortfall estimator
#' \eqn{\widehat{\theta}_T^W}.
#'
#' The weighted fitting procedure computes the sandwich variance quantity
#'
#' \deqn{
#' (\widehat{\Sigma}_T^W)^{-1}
#' \widehat{\Omega}_T^W
#' (\widehat{\Sigma}_T^W)^{-1}.
#' }
#'
#' The estimated covariance matrix of
#' \eqn{\widehat{\theta}_T^W} is the sandwich variance quantity divided by
#' the sample size \eqn{T}. Standard errors are obtained from the square roots
#' of the diagonal entries of this covariance matrix.
#'
#' For each weighted expected shortfall coefficient, the function reports the
#' coefficient estimate, standard error, two-sided normal-approximation
#' confidence interval, and two-sided p-value for testing whether the
#' coefficient is equal to zero.
#'
#' The confidence level is controlled by \code{level}. The variance estimator
#' used during model fitting is determined by \code{covMethod} and
#' \code{truncationMethod} in \code{\link{estools.ts}}.
#'
#' If \code{covariance = TRUE}, the estimated coefficient covariance matrix
#' and the \code{Sigma_hat} and \code{Omega_hat} components from the weighted
#' fitting procedure are also returned.
#'
#' @return
#' An object of class \code{"summary.estools.ts"}. It is a list containing:
#' \describe{
#'   \item{call}{
#'     matched call used to fit the original \code{estools.ts} object.
#'   }
#'   \item{coefficients}{
#'     matrix containing the weighted expected shortfall coefficient estimates,
#'     standard errors, lower and upper confidence limits, and p-values.
#'     The columns are \code{Estimate}, \code{Std. Error}, \code{Lower CI},
#'     \code{Upper CI}, and \code{P-value}.
#'   }
#'   \item{beta.hat}{
#'     estimated coefficients of the conditional quantile model.
#'   }
#'   \item{tau}{
#'     quantile level used for estimation.
#'   }
#'   \item{qmodel}{
#'     autoregressive model used for the conditional quantile.
#'   }
#'   \item{emodel}{
#'     autoregressive model used for the conditional expected shortfall.
#'   }
#'   \item{covMethod}{
#'     conditional variance estimation method used during model fitting.
#'   }
#'   \item{truncationMethod}{
#'     truncated-second-moment estimation method used during model fitting.
#'   }
#'   \item{level}{
#'     confidence level used for the reported confidence intervals.
#'   }
#'   \item{nobs}{
#'     number of observations used for estimation.
#'   }
#'   \item{cov}{
#'     estimated covariance matrix of the weighted expected shortfall
#'     coefficients. Included only when \code{covariance = TRUE}.
#'   }
#'   \item{Sigma_hat}{
#'     \code{Sigma} matrix returned by the weighted fitting procedure.
#'     Included only when \code{covariance = TRUE}.
#'   }
#'   \item{Omega_hat}{
#'     \code{Omega} matrix returned by the weighted fitting procedure.
#'     Included only when \code{covariance = TRUE}.
#'   }
#' }
#'
#' @references
#' Cai, P., Tan, K. M., and He, X.
#' "Estimation and Inference for the Joint Autoregressive
#' Quantile-Expected Shortfall Model."
#'
#' @seealso
#' \code{\link{estools.ts}}
#'
#' @examples
#' set.seed(2026)
#' y <- as.numeric(
#'   stats::arima.sim(
#'     model = list(ar = 0.3),
#'     n = 300
#'   )
#' )
#'
#' fit <- estools.ts(
#'   y,
#'   tau = 0.1,
#'   qmodel = "SAV",
#'   emodel = "SAV",
#'   REP = 1,
#'   MAXITER = 200,
#'   ini_spec = TRUE,
#'   order = 3
#' )
#'
#' summary(fit)
#'
#' summary(
#'   fit,
#'   level = 0.90,
#'   covariance = TRUE
#' )
#'
#' @export
summary.estools.ts <- function(object,
                              level = 0.95,
                              covariance = FALSE,
                              ...){

  if(!inherits(object, "estools.ts")){
    stop("'object' must be a estools.ts object")
  }
  if(!is.numeric(level) || length(level) != 1L || is.na(level) || level <= 0 || level >= 1){
    stop("'level' must be a single number between 0 and 1")
  }

  if(!is.logical(covariance) || length(covariance) != 1L || is.na(covariance)){
    stop("'covariance' must be TRUE or FALSE")
  }
  theta_hat <- object$coefficients
  VC <- object$VC
  n <- length(object$quantile.fitted)
  robust_se <- sqrt(diag(VC) / n)
  z_alpha <- stats::qnorm(1 - (1 - level) / 2)
  lower <- theta_hat - z_alpha * robust_se
  upper <- theta_hat + z_alpha * robust_se
  pivot <- theta_hat / robust_se
  p_values <- 2 * stats::pnorm(abs(pivot), lower.tail = FALSE)
  coef_names <- names(theta_hat)
  if(is.null(coef_names)){
    coef_names <- paste0("theta", seq_along(theta_hat) - 1L)
  }
  coefficients <- cbind(
    Estimate = theta_hat,
    `Std. Error` = robust_se,
    `Lower CI` = lower,
    `Upper CI` = upper,
    `P-value` = p_values
  )
  rownames(coefficients) <- coef_names
  ans <- list(
    call = object$call,
    coefficients = coefficients,
    beta.hat = object$beta.hat,
    tau = object$tau,
    qmodel = object$qmodel,
    emodel = object$emodel,
    covMethod = object$covMethod,
    truncationMethod = object$truncationMethod,
    level = level,
    nobs = n
  )

  if(covariance){
    ans$cov <- VC / n
    ans$Sigma_hat <- object$weighted.fit$Sigma
    ans$Omega_hat <- object$weighted.fit$Omega
  }

  class(ans) <- "summary.estools.ts"
  ans
}
