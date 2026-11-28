#' Fit Linear Expected Shortfall Regression by i-Rock Approach
#'
#' @description
#' Fits a linear conditional expected shortfall (ES) model by projecting
#' initial conditional ES curves using weighted quantile regression.
#' The method is based on the i-Rock framework introduced in Li et al. (2026).
#'
#' @param formula a two-sided formula with one numeric response. Usual
#'   model-matrix terms, including factors and interactions, are supported.
#'   Offsets are not supported.
#' @param data an optional data frame, list, or environment containing variables
#'   in the formula.
#' @param tau one quantile probability strictly between zero and one.
#' @param tail which tail to model. The default is the lower tail, with
#'   probability mass \code{tau}; the upper tail has mass \code{1 - tau}.
#' @param na.action a missing-data action passed to \code{model.frame()}.
#' @param contrasts optional list of contrasts passed to \code{model.matrix()}.
#' @param delta width of the truncated probability window, in \eqn{(0,1]}.
#' @param ds optional positive target spacing for the probability grid.
#'   \code{NULL} uses the sample-size rule described in Appendix G.1 of Li et al. (2026).
#' @param wins fraction of the left side of the paper's upper-tail probability
#'   window to estimate directly, in \eqn{(0,1]}. The omitted part is filled
#'   with a constant. One disables this extension. Ignored for all-discrete
#'   disjoint bins. For lower-tail fits, this operates after response reflection.
#' @param bins optional positive integer number of marginal bins per continuous
#'   covariate. \code{NULL} uses the simulation rule in Section 5 of Li et al. (2026).
#'   Only used with disjoint bins.
#' @param discrete optional character vector naming additional numeric columns
#'   of the predictor model frame to treat as discrete. Factors, characters,
#'   logicals, and numeric columns with fewer than 13 distinct values are
#'   automatically treated as discrete.
#' @param quantile_method initial global quantile fit for continuous or mixed
#'   designs: \code{"linear"} uses the formula design and \code{"spline"}
#'   adds piecewise linear spline terms for the continuous model-frame columns.
#'   All-discrete disjoint designs instead use group-specific sample quantiles.
#' @param ll_weight logical; use local design weights in the final fit for
#'   continuous/mixed bins or neighborhoods. All-discrete disjoint fits always
#'   use group sizes.
#' @param binning use the paper's disjoint bins (default) or the
#'   overlapping nearest-neighbor extension (\code{"knn"}).
#' @param neighbors number of neighbors for the optional nearest-neighbor fit.
#'   Defaults to the smaller of 500 and the number of retained rows.
#' @param centers number of representative observations for the optional
#'   nearest-neighbor fit. Defaults to the smaller of 5 and the number of rows.
#' @param medoids logical; use CLARANS medoids for the nearest-neighbor fit
#'   instead of sampling representative observations without replacement.
#'
#' @details
#' For continuous or mixed covariates, the method fits local linear regressions
#' of the scores at a grid of probability levels and evaluates them at bin
#' representatives. The final linear ES coefficients minimize the integrated
#' quantile loss of these initial ES curves. With only discrete covariates,
#' initial ES curves are empirical tail means within exact covariate groups.
#' The detailed algorithms are included in Li et al. (2026).
#'
#' @return An object of class \code{"irock"}, with components:
#' \describe{
#'   \item{\code{call}}{the matched call.}
#'   \item{\code{formula}}{the expanded formula.}
#'   \item{\code{terms}}{the fitted terms.}
#'   \item{\code{model}}{the model frame.}
#'   \item{\code{na.action}}{the missing-data action.}
#'   \item{\code{contrasts}}{the predictor contrasts.}
#'   \item{\code{xlevels}}{the factor levels.}
#'   \item{\code{tau}}{the requested probability.}
#'   \item{\code{tail}}{the requested tail.}
#'   \item{\code{coefficients}}{the named final ES coefficients.}
#'   \item{\code{quantile}}{the conditional quantile estimates.}
#'   \item{\code{fitted.values}}{the fitted ES values.}
#'   \item{\code{residuals}}{the ES pseudo-response minus fitted ES.}
#' }
#'
#' The \code{quantile}, \code{fitted.values}, and \code{residuals} vectors
#' are computed at the requested level. Like the stored vectors in an
#' \code{lm} fit, these are in retained model-frame order.
#'
#' @references
#' Li, Y., Zhang, S., and He, X. (2026). Expected Shortfall Regression
#' via Optimization. \emph{Journal of the American Statistical Association}, in press.
#'
#' @seealso
#' \code{\link{summary.irock}} for random weighting bootstrap
#' inference on the expected shortfall coefficients.
#'
#' @examples
#' set.seed(21)
#' dat <- data.frame(x = runif(1000))
#' dat$y <- 2 + dat$x + rnorm(1000)
#'
#' fit <- irock(y ~ x, data = dat, tau = 0.1, bins = 3, ds = 0.02)
#' fit$coefficients
#'
#' @export
irock <- function(formula, data, tau = 0.1, tail = c("lower", "upper"),
                  na.action, contrasts = NULL, delta = 0.5, ds = NULL,
                  wins = 0.5, bins = NULL, discrete = NULL,
                  quantile_method = c("linear", "spline"), ll_weight = TRUE,
                  binning = c("disjoint", "knn"), neighbors = NULL,
                  centers = NULL, medoids = TRUE) {
  cl <- match.call()
  if (missing(formula) || !inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula.", call. = FALSE)
  }
  .irock_number(tau, "tau", upper = 1, closed = FALSE)
  .irock_number(delta, "delta", upper = 1)
  .irock_number(wins, "wins", upper = 1)
  if (!is.null(ds)) .irock_number(ds, "ds", upper = 1)
  if (!is.null(bins)) .irock_count(bins, "bins")
  if (!is.null(neighbors)) .irock_count(neighbors, "neighbors")
  if (!is.null(centers)) .irock_count(centers, "centers")
  .irock_flag(ll_weight, "ll_weight")
  .irock_flag(medoids, "medoids")
  tail <- match.arg(tail)
  binning <- match.arg(binning)
  quantile_method <- match.arg(quantile_method)

  mf_call <- cl[c(1L, match(c("formula", "data", "na.action"),
                            names(cl), nomatch = 0L))]
  mf_call[[1L]] <- quote(stats::model.frame)
  mf_call$drop.unused.levels <- TRUE
  mf <- eval(mf_call, parent.frame())
  tr <- attr(mf, "terms")
  y <- stats::model.response(mf)
  if (!is.numeric(y) || !is.null(dim(y)) || !length(y) ||
      any(!is.finite(y))) {
    stop("The retained response must be a nonempty finite numeric vector.",
         call. = FALSE)
  }
  if (length(attr(tr, "offset")) || !is.null(stats::model.offset(mf))) {
    stop("Offsets are not supported by irock().", call. = FALSE)
  }
  x <- stats::model.matrix(tr, mf, contrasts.arg = contrasts)
  if (!ncol(x) || any(!is.finite(x))) {
    stop("The model matrix must have at least one column and finite entries.",
         call. = FALSE)
  }
  if (nrow(x) <= ncol(x) || !.irock_full_rank(x)) {
    stop("The model matrix must have full column rank and more rows than columns.",
         call. = FALSE)
  }
  predictors <- mf[-attr(tr, "response")]
  parts <- .irock_predictors(predictors, discrete, nrow(x))
  xlevels <- lapply(predictors, function(z) {
    if (is.factor(z)) levels(z) else if (is.character(z)) levels(factor(z))
  })
  xlevels <- Filter(Negate(is.null), xlevels)

  sign_y <- if (tail == "lower") -1 else 1
  # The internal algorithm takes a quantile level on its upper-tail scale.
  alpha <- 1 - tau
  if (alpha <= 0 || alpha >= 1) {
    stop("'tau' is too close to a boundary for the working tail.", call. = FALSE)
  }
  groups <- if (binning == "disjoint") .irock_disjoint(parts, bins) else
    .irock_knn(x, neighbors, centers, medoids)
  ans <- .irock_fit(x, sign_y * y, alpha, parts, delta, ds, wins, bins,
                    quantile_method, ll_weight, binning, neighbors, centers,
                    medoids, groups = groups)
  b <- stats::setNames(sign_y * ans$coefficients, colnames(x))
  q <- sign_y * ans$quantile
  es <- drop(x %*% b)
  score <- .irock_score(y, q, tau, tail == "lower")
  if (any(!is.finite(es)) || any(!is.finite(score - es))) {
    stop("Fitted values or ES residuals are nonfinite; rescale the response.",
         call. = FALSE)
  }
  rn <- row.names(mf)
  structure(list(
    call = cl, formula = stats::formula(tr), terms = tr, model = mf,
    na.action = attr(mf, "na.action"), contrasts = attr(x, "contrasts"),
    xlevels = xlevels, tau = tau, tail = tail, coefficients = b,
    quantile = stats::setNames(q, rn), fitted.values = stats::setNames(es, rn),
    residuals = stats::setNames(score - es, rn)
  ), class = "irock", .irock_control = list(
    delta = delta, ds = ds, wins = wins, bins = bins, discrete = discrete,
    quantile_method = quantile_method, ll_weight = ll_weight,
    binning = binning, neighbors = neighbors, centers = centers,
    medoids = medoids, alpha = alpha,
    representative = if (binning == "knn") groups$representative else NULL
  ))
}
