# Internal estimation on the upper-tail scale of Li, Zhang and He (2026).

.irock_score <- function(y, q, tail_mass, lower = FALSE) {
  if (lower) q + pmin(y - q, 0) / tail_mass else q + pmax(y - q, 0) / tail_mass
}

.irock_grid <- function(n, tau, delta, ds, wins, empirical) {
  a <- (1 - delta) * tau
  b <- tau + delta * (1 - tau)
  if (is.null(ds)) {
    ds <- delta / ceiling(sqrt(70 * n * log(n)))
  } else if (!empirical) {
    ds <- max(ds, 3 / (n * log(n)))
  }
  if (delta / ds > 1e5) {
    stop("The probability grid is too large; increase 'ds'.", call. = FALSE)
  }
  intervals <- function(width) max(1, ceiling(width / ds - 1e-10))
  if (empirical) {
    return(list(level = seq(a, b, length.out = 1L + intervals(b - a)),
                weight = NULL))
  }
  start <- (1 - wins * delta) * tau
  # Midpoint quadrature with a boundary at the target probability. This
  # retains exactly the intended probability mass on either side of tau.
  a_edges <- seq(start, tau, length.out = 1L + intervals(tau - start))
  b_edges <- seq(tau, b, length.out = 1L + intervals(b - tau))
  edges <- c(a_edges, b_edges[-1L])
  widths <- diff(edges)
  level <- edges[-length(edges)] + widths / 2
  if (any(widths <= 0) || any(level <= 0 | level >= 1) ||
      !all(level[seq_len(length(a_edges) - 1L)] < tau) ||
      !all(level[seq.int(length(a_edges), length(level))] > tau)) {
    stop("Probability grid collapsed; adjust 'tau', 'delta' or 'wins'.",
         call. = FALSE)
  }
  # Fold all omitted left cells into the first retained response. Check
  # loss is linear in their integration weights, so no copying is needed.
  widths[1L] <- widths[1L] + (start - a)
  list(level = level, weight = widths)
}

.irock_sorted_quantile <- function(ys, levels, w = NULL) {
  n <- length(ys)
  if (is.null(w) || n == 1L || all(w == w[1L])) {
    return(stats::quantile(ys, levels, type = 6L, names = FALSE))
  } else {
    # Weighted plotting positions reduce to j/(n+1) at unit weights.
    # Use convex interpolation to match the original type-6 rule.
    step <- mean(w)
    positions <- (cumsum(w) - w / 2 + step / 2) / (sum(w) + step)
    stats::approx(positions, ys, xout = levels, rule = 2, ties = "ordered")$y
  }
}

.irock_tail_mean <- function(y, levels, weights = NULL) {
  order_y <- order(y)
  ys <- y[order_y]
  n <- length(ys)
  w <- if (is.null(weights)) rep.int(1, n) else weights[order_y]
  q <- .irock_sorted_quantile(ys, levels, w)
  first <- findInterval(q, ys, left.open = TRUE) + 1L
  magnitude <- max(abs(ys))
  if (magnitude == 0) return(rep.int(0, length(levels)))
  total <- rev(cumsum(rev(w * (ys / magnitude))))
  (total[first] / rev(cumsum(rev(w)))[first]) * magnitude
}

.irock_quantile_design <- function(x, xc, method) {
  if (method == "linear") return(x)
  base <- .irock_unit_design(x)
  for (j in seq_len(ncol(xc))) {
    z <- xc[, j]
    limits <- range(z)
    knots <- unique(stats::quantile(z, c(1 / 3, 2 / 3), names = FALSE))
    knots <- knots[knots > limits[1L] & knots < limits[2L]]
    if (!length(knots)) next
    # Only retain independent additions; the original formula columns are
    # always retained (including its chosen intercept convention).
    extra <- splines::bs(z, knots = knots, Boundary.knots = limits, degree = 1L)
    for (k in seq_len(ncol(extra))) {
      col <- extra[, k, drop = FALSE]
      candidate <- cbind(base, .irock_unit_design(col))
      if (.irock_full_rank(candidate) && nrow(candidate) > ncol(candidate)) {
        x <- cbind(x, col)
        base <- candidate
      }
    }
  }
  x
}

.irock_qr <- function(x, y, tau, method, weights = NULL) {
  if (!is.null(weights)) {
    x <- x * weights
    y <- y * weights
  }
  # A zero design row contributes a constant loss; remove it so the solver
  # does not treat its fixed zero prediction as an additional constraint.
  usable <- rowSums(x != 0) > 0L
  if (!all(usable)) {
    x <- x[usable, , drop = FALSE]
    y <- y[usable]
  }
  # Rescale columns for the solver and then return original-scale coefficients.
  scale <- apply(abs(x), 2L, max)
  scale[scale == 0] <- 1
  x <- sweep(x, 2L, scale, "/")
  b <- quantreg::rq.fit(x, y, tau = tau, method = method)$coefficients / scale
  if (length(b) != ncol(x) || any(!is.finite(b))) {
    stop("Quantile regression returned invalid coefficients.", call. = FALSE)
  }
  b
}

.irock_local <- function(xc, idx, representative, weights = NULL) {
  w <- if (is.null(weights)) rep.int(1, length(idx)) else weights[idx]
  if (!ncol(xc)) return(w / sum(w))
  local <- cbind(1, sweep(xc[idx, , drop = FALSE], 2L, xc[representative, ], "-"))
  if (any(!is.finite(local))) {
    stop("Centered covariates overflowed; rescale the predictors.", call. = FALSE)
  }
  # The first coefficient is identifiable even with aliased slopes: the
  # representative has a row with all centered slope coordinates zero.
  # SVD avoids normal equations and supports the Moore-Penrose solution.
  root_w <- sqrt(w)
  s <- svd(.irock_unit_design(local) * root_w)
  keep <- s$d > max(dim(local)) * max(s$d) * .Machine$double.eps
  root_w * drop(s$u[, keep, drop = FALSE] %*% (s$v[1L, keep] / s$d[keep]))
}

# Here tau is the working quantile level, not the public tail mass.
.irock_fit <- function(x, y, tau, parts, delta, ds, wins, bins,
                       quantile_method, ll_weight, binning, neighbors, centers,
                       medoids, groups = NULL, obs_weights = NULL,
                       return_quantile = TRUE) {
  if (is.null(groups)) {
    groups <- if (binning == "disjoint") .irock_disjoint(parts, bins) else
      .irock_knn(x, neighbors, centers, medoids)
  }
  xr <- x[groups$representative, , drop = FALSE]
  if (!.irock_full_rank(xr)) {
    stop("Bin representatives do not identify all ES coefficients; use more bins or centers, or simplify the formula.",
         call. = FALSE)
  }
  n <- length(y)
  grid <- .irock_grid(n, tau, delta, ds, wins, groups$empirical)
  levels <- grid$level
  m <- nrow(xr)
  J <- length(levels)
  if (as.double(m) * J * (1 + ncol(xr)) > 2e7) {
    stop("The stacked fit would be too large; increase 'ds' or reduce the number of bins or centers.",
         call. = FALSE)
  }
  es <- matrix(NA_real_, m, J)
  quantile <- numeric(n)
  if (groups$empirical) {
    weight <- lengths(groups$index)
    for (j in seq_len(m)) {
      idx <- groups$index[[j]]
      es[j, ] <- .irock_tail_mean(y[idx], levels, obs_weights[idx])
      if (!is.null(obs_weights)) weight[j] <- sum(obs_weights[idx])
      if (return_quantile) {
        sorted <- idx[order(y[idx])]
        quantile[idx] <- .irock_sorted_quantile(y[sorted], tau,
                                                obs_weights[sorted])
      }
    }
  } else {
    xq <- .irock_quantile_design(x, parts$continuous, quantile_method)
    local <- Map(function(idx, rep) {
      .irock_local(groups$continuous, idx, rep, obs_weights)
    },
    groups$index, groups$representative)
    # The reciprocal squared norm of the intercept influence vector equals
    # the unnormalized Schur complement (D.2) for full-rank bins.
    weight <- if (ll_weight) {
      if (is.null(obs_weights)) {
        vapply(local, function(a) 1 / sum(a^2), numeric(1L))
      } else {
        vapply(seq_len(m), function(j) {
          1 / sum(local[[j]]^2 / obs_weights[groups$index[[j]]])
        }, numeric(1L))
      }
    } else if (is.null(obs_weights)) {
      rep.int(1, m)
    } else {
      # Perturb unit bin weights by the average multiplier in each bin.
      vapply(groups$index, function(i) mean(obs_weights[i]), numeric(1L))
    }
    # Process one probability at a time, avoiding an n by J quantile matrix.
    for (t in seq_len(J)) {
      alpha <- levels[t]
      q <- drop(xq %*% .irock_qr(xq, y, alpha, "br", obs_weights))
      z <- .irock_score(y, q, 1 - alpha)
      es[, t] <- vapply(seq_len(m), function(j) {
        sum(local[[j]] * z[groups$index[[j]]])
      }, numeric(1L))
    }
    if (return_quantile) {
      quantile <- drop(xq %*% .irock_qr(xq, y, tau, "br", obs_weights))
    }
  }
  weight <- weight / max(weight)
  w <- if (is.null(grid$weight)) rep(weight, J) else
    as.vector(outer(weight, grid$weight / max(grid$weight)))
  # R stacks ES columns in probability order. Repeat representative rows
  # in the same order. Check loss is positively homogeneous in the weights.
  if (any(!is.finite(es)) || any(!is.finite(quantile)) ||
      any(!is.finite(w) | w <= 0)) {
    stop("Initial estimates or integration weights are invalid; rescale the data or adjust the grid.",
         call. = FALSE)
  }
  X <- xr[rep(seq_len(m), J), , drop = FALSE]
  coefficients <- .irock_qr(X * w, as.vector(es) * w, tau,
                            if (groups$empirical) "br" else "fn")
  list(coefficients = coefficients, quantile = quantile)
}
