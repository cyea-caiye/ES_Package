#' Fit Autoregressive Quantile-Expected Shortfall Models
#'
#' @description
#' Fits a joint autoregressive quantile and expected shortfall model using
#' the two-step estimation procedure and the weighted estimator.
#'
#' @param y A numeric time series.
#' @param tau Quantile level. Must be strictly between 0 and 1.
#' @param qmodel Autoregressive model for the conditional quantile.
#'   One of \code{"SAV"}, \code{"AS"}, \code{"GARCH"},
#'   \code{"SAV12"}, or \code{"AS12"}.
#' @param emodel Autoregressive model for the conditional expected shortfall.
#'   By default, the same model as \code{qmodel}.
#' @param REP Number of repetitions used in the optimization procedure.
#' @param MAXITER Maximum number of iterations for each optimization step.
#' @param covMethod Method used for variance-covariance estimation.
#' @param truncationMethod Method used for truncated conditional variance
#'   estimation.
#'
#' @return An object of class \code{"lesreg.ts"} containing the fitted
#'   quantile and expected shortfall models.
#'
#' @export
lesreg.ts <- function(y,
                      tau = 0.05,
                      qmodel = "SAV",
                      emodel = qmodel,
                      REP = 5,
                      MAXITER = 1000,
                      covMethod = 1,
                      truncationMethod = 1){
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
  REP <- as.integer(REP)
  MAXITER <- as.integer(MAXITER)
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
    order = NULL,
    ini_spec = FALSE
  )

  weights <- 1 / (fit.unweighted$denom^2)

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
    warmstart = NULL,
    covMethod = covMethod,
    truncationMethod = truncationMethod,
    q_hat = fit.unweighted$q_hat
  )

  fit <- list(
    coefficients = fit.weighted$res[, "c_est"],
    beta.hat = fit.unweighted$beta_hat,
    theta.hat = fit.weighted$res[, "c_est"],
    theta.unweighted = fit.unweighted$res[, "c_est"],
    quantile.fitted = fit.unweighted$q_hat,
    es.fitted = fit.weighted$ES,
    fitted.values = fit.weighted$ES,
    residuals = fit.weighted$residuals,
    response.residuals = y - fit.weighted$ES,
    pseudo.response = fit.unweighted$pseudo_y,
    weights = weights,
    VC = fit.weighted$VC,
    ci = fit.weighted$res[ , c("ci_low", "ci_up"), drop = FALSE],
    tau = tau,
    qmodel = qmodel,
    emodel = emodel,
    REP = REP,
    MAXITER = MAXITER,
    covMethod = covMethod,
    truncationMethod = truncationMethod,
    unweighted.fit = fit.unweighted,
    weighted.fit = fit.weighted,
    call = match.call()
  )

  class(fit) <- "lesreg.ts"
  fit
}
