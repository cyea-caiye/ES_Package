#' Summarizing Autoregressive Expected Shortfall Regression Fits
#'
#' @description
#' Summarizes a fitted autoregressive expected shortfall (ES) regression model
#' and provides standard errors, confidence intervals, and \eqn{p}-values for
#' the final weighted ES regression coefficients.
#'
#' @param x an object of class \code{"summary.esreg.ts"}, usually a result of a
#'   call to \code{summary.esreg.ts}.
#' @param object an object of class \code{"esreg.ts"}, usually a result of a
#'   call to \code{\link{esreg.ts}}.
#' @param level confidence level for the reported confidence intervals. This
#'   must be a single number strictly between 0 and 1.
#' @param covariance logical flag indicating whether the estimated covariance
#'   matrix and its sandwich components should be returned.
#' @param digits the number of significant digits to use when printing.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' A fitted \code{"esreg.ts"} object contains the first-step conditional
#' quantile regression coefficients together with the unweighted and weighted
#' expected shortfall regression coefficients. The summary method performs
#' statistical inference for the final weighted expected shortfall regression
#' coefficients. See \code{\link{esreg.ts}} for details on the underlying
#' time-series two-step estimation procedure.
#'
#' Standard errors are computed using the sandwich covariance matrix obtained
#' during model fitting. The construction of this matrix depends on the
#' \code{covMethod} and \code{truncationMethod} specified in the call to
#' \code{\link{esreg.ts}}.
#'
#' The reported \eqn{p}-values are computed using a standard normal
#' approximation.
#'
#' @return
#' \code{summary.esreg.ts} returns an object of class
#' \code{"summary.esreg.ts"} containing summary statistics and inference
#' results for the fitted autoregressive expected shortfall regression model.
#' The documented components are:
#' \describe{
#'   \item{\code{coefficients}}{a matrix containing the estimated weighted
#'   expected shortfall regression coefficients, their standard errors, lower
#'   and upper confidence limits, and \eqn{p}-values.}
#'   \item{\code{beta.hat}}{the estimated coefficients from the first-step
#'   conditional quantile regression fit.}
#'   \item{\code{tau}}{the quantile level used for estimation.}
#'   \item{\code{qmodel}}{the autoregressive specification used for the
#'   conditional quantile.}
#'   \item{\code{emodel}}{the autoregressive specification used for the
#'   conditional expected shortfall.}
#'   \item{\code{covMethod}}{the conditional variance estimation method used
#'   in model fitting.}
#'   \item{\code{truncationMethod}}{the truncated second-moment estimation
#'   method used in model fitting.}
#'   \item{\code{level}}{the confidence level used for the reported confidence
#'   intervals.}
#'   \item{\code{nobs}}{the number of observations used in the fitted model.}
#'   \item{\code{cov}}{the estimated covariance matrix of the weighted
#'   expected shortfall regression coefficients, if
#'   \code{covariance = TRUE}.}
#'   \item{\code{Sigma_hat}}{the estimated outer matrix in the sandwich
#'   covariance calculation, if \code{covariance = TRUE}.}
#'   \item{\code{Omega_hat}}{the estimated middle matrix in the sandwich
#'   covariance calculation, if \code{covariance = TRUE}.}
#' }
#'
#' @seealso
#' \code{\link{esreg.ts}}
#'
#' @examples
#' set.seed(2026)
#' T <- 5000
#' tau <- 0.1
#' alpha <- c(0.45, 0.65, -0.25, -0.35, 0.8, 0.35, 0.25)
#'
#' y <- numeric(T)
#' f <- numeric(T)
#' f[1] <- alpha[1] / (1 - alpha[2])
#' y[1] <- f[1] + alpha[5] * rnorm(1)
#'
#' for(t in 2:T){
#'   y.pos <- y[t - 1] * (y[t - 1] >= 0)
#'   y.neg <- y[t - 1] * (y[t - 1] < 0)
#'   f[t] <- alpha[1] + alpha[2] * f[t - 1] +
#'     alpha[3] * y.pos - alpha[4] * y.neg
#'   scale.t <- alpha[5] + alpha[6] * y.pos - alpha[7] * y.neg
#'   y[t] <- f[t] + scale.t * rnorm(1)
#' }
#'
#' fit <- esreg.ts(
#'   y = y,
#'   tau = tau,
#'   qmodel = "AS12",
#'   emodel = "AS12",
#'   REP = 5,
#'   MAXITER = 1000,
#'   ini_spec = FALSE
#' )
#'
#' ## Summarize the fitted model
#' summary(fit)
#'
#' ## Return the estimated covariance matrix and sandwich components
#' summary(fit, covariance = TRUE)
#'
#' @export
summary.esreg.ts <- function(object,
                              level = 0.95,
                              covariance = FALSE,
                              ...){

  if(!inherits(object, "esreg.ts")){
    stop("'object' must be a esreg.ts object")
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

  class(ans) <- "summary.esreg.ts"
  ans
}
