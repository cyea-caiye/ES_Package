#' Autoregressive Expected Shortfall Regression
#'
#' @description
#' Fits an autoregressive expected shortfall (ES) regression model for a
#' univariate time series using a two-step regression approach. The function
#' also fits the corresponding autoregressive quantile regression model, which
#' is used in the first step of estimation.
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
#' @param y a numeric vector representing a univariate time series.
#' @param tau the quantile level associated with the ES to be estimated. This
#'   must be a single number strictly between 0 and 1.
#' @param qmodel character string specifying the autoregressive model for the
#'   conditional quantile. Available choices are \code{"SAV"}, \code{"AS"},
#'   \code{"GARCH"}, \code{"SAV12"}, and \code{"AS12"}.
#' @param emodel character string specifying the autoregressive model for the
#'   conditional expected shortfall. Available choices are \code{"SAV"},
#'   \code{"AS"}, \code{"GARCH"}, \code{"SAV12"}, and \code{"AS12"}.
#'   By default, \code{emodel = qmodel}. If \code{ini_spec = TRUE},
#'   \code{qmodel} and \code{emodel} must be the same.
#' @param REP positive integer specifying the number of repetitions used in
#'   the optimization procedure.
#' @param MAXITER positive integer specifying the maximum number of iterations
#'   for each optimization step.
#' @param covMethod integer specifying the method used for conditional
#'   variance estimation. \code{1} uses squared residuals and \code{2} uses
#'   an estimator based on the truncated second moment.
#' @param truncationMethod integer specifying the method used to estimate the
#'   truncated second moment. \code{1} uses the empirical truncated second
#'   moment and \code{2} uses a kernel-density-based numerical approximation.
#' @param ini_spec logical. If \code{FALSE}, random initial values are used in
#'   numerical optimization. If \code{TRUE}, sieve-based initial values are
#'   used, and \code{qmodel} and \code{emodel} must specify the same model.
#' @param order positive integer specifying the sieve order when
#'   \code{ini_spec = TRUE}. If \code{NULL}, it is set to
#'   \code{round(log(length(y)))}.
#'
#' @details
#' Let \eqn{e_t(\boldsymbol{\theta})} denote the conditional expected
#' shortfall at time \eqn{t}. The argument \code{emodel} selects its
#' autoregressive specification from the following models.
#'
#' Define
#' \eqn{Y_t^+ = Y_t\mathbf{1}(Y_t \geq 0)} and
#' \eqn{Y_t^- = Y_t\mathbf{1}(Y_t < 0)}.
#'
#' For \code{"SAV"}, the symmetric absolute value specification is
#' \deqn{
#' e_t(\boldsymbol{\theta})
#' =
#' \theta_0 + \theta_1 e_{t-1}(\boldsymbol{\theta})
#' + \theta_2 |Y_{t-1}|.
#' }
#'
#' For \code{"AS"}, the asymmetric slope specification is
#' \deqn{
#' e_t(\boldsymbol{\theta})
#' =
#' \theta_0 + \theta_1 e_{t-1}(\boldsymbol{\theta})
#' + \theta_2 Y_{t-1}^+ - \theta_3 Y_{t-1}^-.
#' }
#'
#' For \code{"GARCH"}, the specification is
#' \deqn{
#' e_t(\boldsymbol{\theta})
#' =
#' -\sqrt{
#' \theta_0 + \theta_1 e_{t-1}^2(\boldsymbol{\theta})
#' + \theta_2 Y_{t-1}^2
#' }.
#' }
#'
#' For \code{"SAV12"}, the lag-two symmetric absolute value specification is
#' \deqn{
#' e_t(\boldsymbol{\theta})
#' =
#' \theta_0 + \theta_1 e_{t-1}(\boldsymbol{\theta})
#' + \theta_2 |Y_{t-1}| + \theta_3 |Y_{t-2}|.
#' }
#'
#' For \code{"AS12"}, the lag-two asymmetric slope specification is
#' \deqn{
#' e_t(\boldsymbol{\theta})
#' =
#' \theta_0 + \theta_1 e_{t-1}(\boldsymbol{\theta})
#' + \theta_2 Y_{t-1}^+ - \theta_3 Y_{t-1}^-
#' + \theta_4 Y_{t-2}^+ - \theta_5 Y_{t-2}^-.
#' }
#'
#' The conditional quantile \eqn{q_t(\boldsymbol{\beta})} uses the same set
#' of autoregressive specifications, with
#' \eqn{e_t(\boldsymbol{\theta})} and \eqn{\boldsymbol{\theta}} replaced by
#' \eqn{q_t(\boldsymbol{\beta})} and \eqn{\boldsymbol{\beta}},
#' respectively. The quantile and expected shortfall specifications are
#' selected by \code{qmodel} and \code{emodel}.
#'
#' When \code{ini_spec = FALSE}, the conditional quantile and expected
#' shortfall may use different autoregressive specifications. When
#' \code{ini_spec = TRUE}, \code{qmodel} and \code{emodel} must specify the
#' same model.
#'
#' @section Method:
#' The function follows the same two-step expected shortfall regression
#' principle described in \code{\link{lesreg}}. In the time-series setting,
#' the conditional quantile and expected shortfall are recursively defined
#' functions of their parameters and are therefore fitted by numerical
#' optimization.
#'
#' An unweighted autoregressive ES fit is first obtained. Let
#' \eqn{\hat q_t} and \eqn{\hat e_t} denote the fitted conditional quantile
#' and unweighted conditional expected shortfall, respectively. The weights
#' for the final weighted fit are constructed from
#' \deqn{
#' w_t = (\hat q_t - \hat e_t)^{-2}.
#' }
#' The final ES coefficients are then obtained from the corresponding weighted
#' autoregressive ES fit.
#'
#' If \code{ini_spec = FALSE}, random initial values are used in numerical
#' optimization. If \code{ini_spec = TRUE}, sieve-based initial values are
#' used. In this case, the conditional quantile and expected shortfall must
#' have the same autoregressive specification, and the sieve order is
#' controlled by \code{order}.
#'
#' @return
#' \code{esreg.ts} returns an object of class \code{"esreg.ts"}.
#' The functions \code{\link{summary.esreg.ts}} and
#' \code{\link{print.esreg.ts}} can be used to summarize and print the fitted
#' model.
#'
#' An object of class \code{"esreg.ts"} is a list containing at least the
#' following components:
#' \describe{
#'   \item{\code{coefficients}}{estimated coefficients from the final weighted
#'   expected shortfall regression fit.}
#'   \item{\code{beta.hat}}{estimated coefficients from the first-step
#'   conditional quantile regression fit.}
#'   \item{\code{theta.hat}}{estimated coefficients from the final weighted
#'   expected shortfall regression fit.}
#'   \item{\code{theta.unweighted}}{estimated coefficients from the unweighted
#'   expected shortfall regression fit.}
#'   \item{\code{quantile.fitted}}{fitted conditional quantile values.}
#'   \item{\code{es.fitted}}{fitted conditional expected shortfall values from
#'   the final weighted fit.}
#'   \item{\code{fitted.values}}{fitted conditional expected shortfall values
#'   from the final weighted fit.}
#'   \item{\code{residuals}}{residuals from the final weighted expected
#'   shortfall regression fit.}
#'   \item{\code{response.residuals}}{the observed response minus the fitted
#'   conditional expected shortfall values.}
#'   \item{\code{pseudo.response}}{the pseudo-response used for expected
#'   shortfall estimation.}
#'   \item{\code{weights}}{weights used in the final weighted expected
#'   shortfall regression fit.}
#'   \item{\code{VC}}{the estimated sandwich matrix used for inference on the
#'   weighted expected shortfall coefficients.}
#'   \item{\code{ci}}{confidence intervals computed during model fitting.}
#'   \item{\code{tau}}{the quantile level used for estimation.}
#'   \item{\code{qmodel}}{the autoregressive specification used for the
#'   conditional quantile.}
#'   \item{\code{emodel}}{the autoregressive specification used for the
#'   conditional expected shortfall.}
#'   \item{\code{REP}}{the number of optimization repetitions used.}
#'   \item{\code{MAXITER}}{the maximum number of iterations used in each
#'   optimization step.}
#'   \item{\code{covMethod}}{the conditional variance estimation method used.}
#'   \item{\code{truncationMethod}}{the truncated second-moment estimation
#'   method used.}
#'   \item{\code{ini_spec}}{logical flag indicating whether sieve-based
#'   initial values were used.}
#'   \item{\code{order}}{the sieve order used when applicable.}
#'   \item{\code{unweighted.fit}}{the internal unweighted expected shortfall
#'   fit.}
#'   \item{\code{weighted.fit}}{the internal weighted expected shortfall fit.}
#'   \item{\code{call}}{the call that produced the fitted object.}
#' }
#'
#' @seealso
#' \code{\link{summary.esreg.ts}}, \code{\link{lesreg}}
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
#' es.eta <- -stats::dnorm(stats::qnorm(tau)) / tau
#' theta.true <- c(
#'   alpha[1] + (1 - alpha[2]) * alpha[5] * es.eta,
#'   alpha[2],
#'   alpha[3] + alpha[6] * es.eta,
#'   alpha[4] + alpha[7] * es.eta,
#'   -alpha[2] * alpha[6] * es.eta,
#'   -alpha[2] * alpha[7] * es.eta
#' )
#'
#' ## Fit the model using random initial values
#' fit <- esreg.ts(
#'   y = y,
#'   tau = tau,
#'   qmodel = "AS12",
#'   emodel = "AS12",
#'   REP = 5,
#'   MAXITER = 1000,
#'   ini_spec = FALSE
#' )
#' theta.true
#' fit$coefficients
#'
#' ## Fit the model using sieve-based initial values
#' fit.ini <- esreg.ts(
#'   y = y,
#'   tau = tau,
#'   qmodel = "AS12",
#'   emodel = "AS12",
#'   REP = 5,
#'   MAXITER = 1000,
#'   ini_spec = TRUE
#' )
#' fit.ini$coefficients
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
