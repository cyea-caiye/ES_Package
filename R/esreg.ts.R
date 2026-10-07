#' Fitting Autoregressive Quantile-Expected Shortfall Models
#'
#' @description
#' \code{esreg.ts} fits autoregressive conditional quantile and expected
#' shortfall models for time-series data using a weighted two-step estimation
#' procedure. The conditional quantile is first estimated by nonlinear quantile
#' regression. An unweighted expected shortfall estimator is then obtained
#' using a constructed pseudo-response, followed by a weighted expected
#' shortfall estimator.
#'
#' @usage
#' esreg.ts(
#'   y,
#'   tau = 0.05,
#'   qmodel = "SAV",
#'   emodel = qmodel,
#'   REP = 5,
#'   MAXITER = 1000,
#'   covMethod = 1,
#'   truncationMethod = 1,
#'   ini_spec = FALSE,
#'   order = NULL
#' )
#'
#' @param y a numeric time series.
#' @param tau quantile level. It must be a single number strictly between
#'   0 and 1.
#' @param qmodel autoregressive model used for the conditional quantile.
#'   Available choices are \code{"SAV"}, \code{"AS"}, \code{"GARCH"},
#'   \code{"SAV12"}, and \code{"AS12"}.
#' @param emodel autoregressive model used for the conditional expected
#'   shortfall. Available choices are \code{"SAV"}, \code{"AS"},
#'   \code{"GARCH"}, \code{"SAV12"}, and \code{"AS12"}. By default,
#'   \code{emodel = qmodel}.
#' @param REP positive integer giving the number of repeated optimization
#'   steps.
#' @param MAXITER positive integer giving the maximum number of iterations
#'   used in each optimization step.
#' @param covMethod integer specifying the conditional variance estimator used
#'   for inference. \code{1} uses squared residuals. \code{2} uses the
#'   conditional variance estimator based on the truncated second moment.
#' @param truncationMethod integer specifying how the truncated second moment
#'   is estimated when \code{covMethod = 2}. It must be either \code{1} or
#'   \code{2}. Method \code{1} uses the empirical truncated second moment,
#'   whereas method \code{2} uses a kernel-density based numerical
#'   integration.
#' @param ini_spec logical. If \code{TRUE}, sieve-based initial values are used
#'   for optimization, and \code{qmodel} and \code{emodel} must specify the
#'   same model. If \code{FALSE}, random initial values are used.
#' @param order positive integer giving the sieve order used when
#'   \code{ini_spec = TRUE}. If \code{NULL}, it is set to
#'   \code{round(log(length(y)))}.
#'
#' @details
#' Let \eqn{q_t(\beta)} and \eqn{e_t(\theta)} denote the conditional quantile
#' and conditional expected shortfall functions, respectively. The current
#' implementation provides the following autoregressive specifications.
#'
#' For \code{"SAV"}, the conditional quantile recursion is
#'
#' \deqn{
#' q_t(\beta)
#' =
#' \beta_0
#' +
#' \beta_1 q_{t-1}(\beta)
#' +
#' \beta_2 |Y_{t-1}|,
#' }
#'
#' and the conditional expected shortfall recursion is
#'
#' \deqn{
#' e_t(\theta)
#' =
#' \theta_0
#' +
#' \theta_1 e_{t-1}(\theta)
#' +
#' \theta_2 |Y_{t-1}|.
#' }
#'
#' For \code{"AS"}, the conditional quantile recursion implemented in the
#' package is
#'
#' \deqn{
#' q_t(\beta)
#' =
#' \beta_0
#' +
#' \beta_1 q_{t-1}(\beta)
#' +
#' \beta_2 Y_{t-1} 1(Y_{t-1}>0)
#' -
#' \beta_3 Y_{t-1} 1(Y_{t-1}<0),
#' }
#'
#' and the conditional expected shortfall recursion is
#'
#' \deqn{
#' e_t(\theta)
#' =
#' \theta_0
#' +
#' \theta_1 e_{t-1}(\theta)
#' +
#' \theta_2 Y_{t-1} 1(Y_{t-1}>0)
#' -
#' \theta_3 Y_{t-1} 1(Y_{t-1}<0).
#' }
#'
#' For \code{"GARCH"}, the conditional quantile recursion is
#'
#' \deqn{
#' q_t(\beta)
#' =
#' -\sqrt{
#' \beta_0
#' +
#' \beta_1 q_{t-1}^2(\beta)
#' +
#' \beta_2 Y_{t-1}^2
#' },
#' }
#'
#' and the conditional expected shortfall recursion is
#'
#' \deqn{
#' e_t(\theta)
#' =
#' -\sqrt{
#' \theta_0
#' +
#' \theta_1 e_{t-1}^2(\theta)
#' +
#' \theta_2 Y_{t-1}^2
#' }.
#' }
#'
#' For \code{"SAV12"}, the conditional quantile recursion implemented in the
#' package is
#'
#' \deqn{
#' q_t(\beta)
#' =
#' \beta_0
#' +
#' \beta_1 q_{t-1}(\beta)
#' +
#' \beta_2 |Y_{t-1}|
#' +
#' \beta_3 |Y_{t-2}|,
#' }
#'
#' and the conditional expected shortfall recursion is
#'
#' \deqn{
#' e_t(\theta)
#' =
#' \theta_0
#' +
#' \theta_1 e_{t-1}(\theta)
#' +
#' \theta_2 |Y_{t-1}|
#' +
#' \theta_3 |Y_{t-2}|.
#' }
#'
#' For \code{"AS12"}, the conditional quantile recursion implemented in the
#' package is
#'
#' \deqn{
#' q_t(\beta)
#' =
#' \beta_0
#' +
#' \beta_1 q_{t-1}(\beta)
#' +
#' \beta_2 Y_{t-1} 1(Y_{t-1}>0)
#' -
#' \beta_3 Y_{t-1} 1(Y_{t-1}<0)
#' +
#' \beta_4 Y_{t-2} 1(Y_{t-2}>0)
#' -
#' \beta_5 Y_{t-2} 1(Y_{t-2}<0),
#' }
#'
#' and the conditional expected shortfall recursion is
#'
#' \deqn{
#' e_t(\theta)
#' =
#' \theta_0
#' +
#' \theta_1 e_{t-1}(\theta)
#' +
#' \theta_2 Y_{t-1} 1(Y_{t-1}>0)
#' -
#' \theta_3 Y_{t-1} 1(Y_{t-1}<0)
#' +
#' \theta_4 Y_{t-2} 1(Y_{t-2}>0)
#' -
#' \theta_5 Y_{t-2} 1(Y_{t-2}<0).
#' }
#'
#' The estimation procedure follows the weighted two-step autoregressive
#' expected shortfall regression procedure. Define
#'
#' \deqn{
#' S_0(\beta,\theta;Y_t,X_t)
#' =
#' e_t(\theta)
#' -
#' q_t(\beta)
#' -
#' \frac{
#' 1\{Y_t \leq q_t(\beta)\}
#' \{Y_t-q_t(\beta)\}
#' }{\tau}.
#' }
#'
#' In the first step, the conditional quantile parameter is estimated from
#'
#' \deqn{
#' \widehat{\beta}_T
#' \in
#' \arg\min_{\beta \in B}
#' \frac{1}{T}
#' \sum_{t=1}^{T}
#' \rho_{\tau}
#' \{Y_t-q_t(\beta)\},
#' }
#'
#' where
#'
#' \deqn{
#' \rho_{\tau}(u)
#' =
#' u[\tau-1\{u\leq0\}].
#' }
#'
#' Given \eqn{\widehat{\beta}_T}, the unweighted expected shortfall estimator
#' is obtained from
#'
#' \deqn{
#' \widehat{\theta}_T
#' \in
#' \arg\min_{\theta \in \Theta}
#' \frac{1}{2T}
#' \sum_{t=1}^{T}
#' S_0^2(
#' \widehat{\beta}_T,
#' \theta;
#' Y_t,
#' X_t).
#' }
#'
#' Equivalently, this is a nonlinear least-squares problem with the
#' pseudo-response
#'
#' \deqn{
#' q_t(\widehat{\beta}_T)
#' +
#' \frac{
#' 1\{Y_t \leq q_t(\widehat{\beta}_T)\}
#' \{Y_t-q_t(\widehat{\beta}_T)\}
#' }{\tau}.
#' }
#'
#' The final weighted expected shortfall estimator is
#'
#' \deqn{
#' \widehat{\theta}_T^W
#' \in
#' \arg\min_{\theta \in \Theta}
#' \frac{1}{2T}
#' \sum_{t=1}^{T}
#' \frac{
#' S_0^2(
#' \widehat{\beta}_T,
#' \theta;
#' Y_t,
#' X_t)
#' }{
#' \{
#' q_t(\widehat{\beta}_T)
#' -
#' e_t(\widehat{\theta}_T)
#' \}^2
#' }.
#' }
#'
#' The implementation constructs the weights from the fitted conditional
#' quantile and the fitted unweighted expected shortfall. Very small values of
#' their difference are stabilized internally before the inverse-square
#' weights are computed.
#'
#' When \code{ini_spec = FALSE}, random initial parameter vectors are used in
#' the optimization. When \code{ini_spec = TRUE}, sieve-based initial values
#' are used, with the sieve order controlled by \code{order}.
#'
#' The numerical optimization uses \code{\link[stats]{optim}} with
#' L-BFGS-B and Nelder-Mead steps.
#'
#' Inference for the final weighted expected shortfall estimator is based on
#' the sandwich variance calculation implemented during the weighted fitting
#' step. See \code{\link{summary.esreg.ts}} for the reported standard
#' errors, confidence intervals, and p-values.
#'
#' @return
#' An object of class \code{"esreg.ts"}. It is a list containing:
#' \describe{
#'   \item{coefficients}{
#'     estimated coefficients of the final weighted expected shortfall model.
#'   }
#'   \item{beta.hat}{
#'     estimated coefficients of the conditional quantile model.
#'   }
#'   \item{theta.unweighted}{
#'     estimated coefficients from the unweighted expected shortfall fit.
#'   }
#'   \item{quantile.fitted}{
#'     fitted conditional quantile values.
#'   }
#'   \item{fitted.values}{
#'     fitted conditional expected shortfall values from the final weighted
#'     fit.
#'   }
#'   \item{residuals}{
#'     residuals from the pseudo-response regression in the final weighted
#'     expected shortfall fit.
#'   }
#'   \item{response.residuals}{
#'     observed responses minus the fitted conditional expected shortfall
#'     values.
#'   }
#'   \item{pseudo.response}{
#'     pseudo-response used for expected shortfall estimation.
#'   }
#'   \item{weights}{
#'     weights used in the final weighted expected shortfall fit.
#'   }
#'   \item{VC}{
#'     sandwich variance quantity returned by the weighted fitting procedure.
#'     The estimated coefficient covariance matrix used by
#'     \code{summary.esreg.ts} is \code{VC / length(y)}.
#'   }
#'   \item{ci}{
#'     95 percent confidence intervals for the final weighted expected
#'     shortfall coefficients computed during model fitting.
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
#'   \item{REP}{
#'     number of repeated optimization steps.
#'   }
#'   \item{MAXITER}{
#'     maximum number of iterations used in each optimization step.
#'   }
#'   \item{covMethod}{
#'     conditional variance estimation method used for inference.
#'   }
#'   \item{truncationMethod}{
#'     truncated-second-moment estimation method used for inference.
#'   }
#'   \item{ini_spec}{
#'     logical flag indicating whether sieve-based initial values were used.
#'   }
#'   \item{order}{
#'     sieve order used when applicable.
#'   }
#'   \item{unweighted.fit}{
#'     detailed output from the unweighted expected shortfall fitting step.
#'   }
#'   \item{weighted.fit}{
#'     detailed output from the weighted expected shortfall fitting step.
#'   }
#'   \item{call}{
#'     matched function call.
#'   }
#' }
#'
#'
#' @seealso
#' \code{\link{summary.esreg.ts}},
#' \code{\link{lesreg}}
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
#' fit <- esreg.ts(
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
#' fit
#' summary(fit)
#'
#' @export
esreg.ts <- function(y,
                      tau = 0.05,
                      qmodel = "SAV",
                      emodel = qmodel,
                      REP = 5,
                      MAXITER = 1000,
                      covMethod = 1,
                      truncationMethod = 1,
                      ini_spec = FALSE,
                      order = NULL){
  if(!is.numeric(y)){
    stop("'y' must be numeric")
  }
  y <- as.numeric(y)
  if(any(!is.finite(y))){
    stop("'y' must contain only finite values")
  }
  if(!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau <= 0 || tau >= 1){
    stop("'tau' must be a single number strictly between 0 and 1")
  }
  if(!is.numeric(REP) || length(REP) != 1L || is.na(REP) || REP < 1){
    stop("'REP' must be a positive integer")
  }
  if(!is.numeric(MAXITER) || length(MAXITER) != 1L || is.na(MAXITER) || MAXITER < 1){
    stop("'MAXITER' must be a positive integer")
  }
  if(!is.numeric(covMethod) || length(covMethod) != 1L || is.na(covMethod) ||
     covMethod != as.integer(covMethod) || !covMethod %in% c(1L, 2L)){
    stop("'covMethod' must be either 1 or 2")
  }
  if(!is.numeric(truncationMethod) || length(truncationMethod) != 1L || is.na(truncationMethod) ||
     truncationMethod != as.integer(truncationMethod) || !truncationMethod %in% c(1L, 2L)){
    stop("'truncationMethod' must be either 1 or 2")
  }
  if(!is.logical(ini_spec) || length(ini_spec) != 1L || is.na(ini_spec)){
    stop("'ini_spec' must be TRUE or FALSE")
  }

  if(!is.null(order)){
    if(!is.numeric(order) || length(order) != 1L || is.na(order) || order < 1 || order != as.integer(order)){
      stop("'order' must be a positive integer")
    }
    order <- as.integer(order)
  }
  if(ini_spec && is.null(order)){
    order <- as.integer(round(log(length(y)), 0))
  }
  REP <- as.integer(REP)
  MAXITER <- as.integer(MAXITER)
  covMethod <- as.integer(covMethod)
  truncationMethod <- as.integer(truncationMethod)
  models <- c(
    SAV   = 1L,
    AS    = 2L,
    GARCH = 3L,
    SAV12 = 4L,
    AS12  = 5L
  )
  if(!is.character(qmodel) || length(qmodel) != 1L || is.na(qmodel)){
    stop("'qmodel' must be a single character string")
  }
  if(!is.character(emodel) || length(emodel) != 1L || is.na(emodel)){
    stop("'emodel' must be a single character string")
  }
  qmodel <- toupper(qmodel)
  emodel <- toupper(emodel)
  if(!qmodel %in% names(models)){
    stop("'qmodel' must be one of: ", paste(names(models), collapse = ", "))
  }
  if(!emodel %in% names(models)){
    stop("'emodel' must be one of: ", paste(names(models), collapse = ", "))
  }
  if(ini_spec && qmodel != emodel){
    warning("When 'ini_spec = TRUE', 'qmodel' and 'emodel' must be the same.", call. = FALSE)
    stop("Please use the same model for 'qmodel' and 'emodel' when 'ini_spec = TRUE'.", call. = FALSE)
  }
  qmodel.id <- unname(models[qmodel])
  emodel.id <- unname(models[emodel])

  fit.unweighted <- ESEst_Inf_fair(
    qmodel = qmodel.id,
    model = emodel.id,
    pval = tau,
    data = y,
    weights = rep(1, length(y)),
    predict = FALSE,
    REP = REP,
    MAXITER = MAXITER,
    covMethod = covMethod,
    truncationMethod = truncationMethod,
    order = order,
    ini_spec = ini_spec
  )

  weights <- 1 / (fit.unweighted$denom^2)
  if(ini_spec){
    warmstart <- fit.unweighted$res[, "c_est"]
  } else{
    warmstart <- NULL
  }

  fit.weighted <- ESEst_Inf_wls_fair(
    pseudo_y = fit.unweighted$pseudo_y,
    emp_es = fit.unweighted$emp_es,
    model = emodel.id,
    pval = tau,
    data = y,
    weights = weights,
    predict = FALSE,
    REP = REP,
    MAXITER = MAXITER,
    warmstart = warmstart,
    covMethod = covMethod,
    truncationMethod = truncationMethod,
    q_hat = fit.unweighted$q_hat
  )

  fit <- list(
    coefficients = fit.weighted$res[, "c_est"],
    beta.hat = fit.unweighted$beta_hat,
    theta.unweighted = fit.unweighted$res[, "c_est"],
    quantile.fitted = fit.unweighted$q_hat,
    fitted.values = fit.weighted$ES,
    residuals = fit.weighted$residuals,
    response.residuals = y - fit.weighted$ES,
    pseudo.response = fit.unweighted$pseudo_y,
    weights = weights,
    VC = fit.weighted$VC,
    ci = fit.weighted$res[, c("ci_low", "ci_up"), drop = FALSE],
    tau = tau,
    qmodel = qmodel,
    emodel = emodel,
    REP = REP,
    MAXITER = MAXITER,
    covMethod = covMethod,
    truncationMethod = truncationMethod,
    ini_spec = ini_spec,
    order = order,
    unweighted.fit = fit.unweighted,
    weighted.fit = fit.weighted,
    call = match.call()
  )

  class(fit) <- "esreg.ts"
  fit
}
