# Internal model validation, covariate groups and bin representatives.

.irock_number <- function(x, name, upper = Inf, closed = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || !is.null(dim(x)) ||
      !is.finite(x) || x <= 0 || x > upper || (!closed && x == upper)) {
    stop(sprintf("'%s' must be one finite number in (0, %s%s.",
                 name, format(upper), if (closed) "]" else ")"),
         call. = FALSE)
  }
}

.irock_count <- function(x, name) {
  .irock_number(x, name, upper = .Machine$integer.max)
  if (x != floor(x)) {
    stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
  }
}

.irock_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || !is.null(dim(x)) || is.na(x)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
}

.irock_unit_design <- function(x) {
  scale <- apply(abs(x), 2L, max)
  scale[scale == 0] <- 1
  sweep(x, 2L, scale, "/")
}

.irock_full_rank <- function(x) {
  nrow(x) >= ncol(x) &&
    qr(.irock_unit_design(x), tol = 1e-10)$rank == ncol(x)
}

.irock_predictors <- function(mf, discrete, n) {
  continuous <- groups <- list()
  if (!is.null(discrete) &&
      (!is.character(discrete) || anyNA(discrete) ||
       !all(discrete %in% names(mf)))) {
    stop("'discrete' must name columns of the predictor model frame.",
         call. = FALSE)
  }
  for (j in seq_along(mf)) {
    z <- mf[[j]]
    name <- names(mf)[j]
    valid <- function(v) is.numeric(v) || is.factor(v) ||
      is.character(v) || is.logical(v)
    if (!valid(z) || (is.array(z) && !is.matrix(z))) {
      stop(sprintf("Unsupported predictor type for '%s'.", name), call. = FALSE)
    }
    # Matrix terms such as poly(x, 2) contribute their model-frame columns.
    for (v in if (is.matrix(z)) as.data.frame(z) else list(z)) {
      if (anyNA(v) || (is.numeric(v) && any(!is.finite(v)))) {
        stop(sprintf("Predictor '%s' has missing or nonfinite values.", name),
             call. = FALSE)
      }
      u <- unique(v)
      if (!is.numeric(v) || name %in% discrete || length(u) < 13L) {
        groups[[length(groups) + 1L]] <- match(v, u)
      } else {
        continuous[[length(continuous) + 1L]] <- v
      }
    }
  }
  list(
    continuous = if (length(continuous)) do.call(cbind, continuous) else
      matrix(numeric(), n, 0L),
    discrete = groups
  )
}

.irock_groups <- function(ids, n) {
  if (!length(ids)) return(list(seq_len(n)))
  if (length(ids) == 1L) {
    key <- ids[[1L]]
  } else {
    key <- do.call(paste, c(unname(ids), sep = ":"))
  }
  # Only form occupied groups; never enumerate a Cartesian product of levels.
  unname(split(seq_len(n), match(key, unique(key))))
}

.irock_disjoint <- function(parts, bins) {
  xc <- parts$continuous
  n <- nrow(xc)
  d <- ncol(xc)
  if (!d) {
    idx <- .irock_groups(parts$discrete, n)
    return(list(index = idx, representative = vapply(idx, "[", 1L, FUN.VALUE = 1L),
                continuous = xc, empirical = TRUE))
  }
  if (is.null(bins)) {
    bins <- min(n, ceiling(1.6 * sqrt(d) * (sqrt(n) / log(n))^(1 / d)))
  }
  if (bins > n) stop("'bins' cannot exceed the number of retained rows.",
                     call. = FALSE)
  cells <- midpoint <- vector("list", d)
  probs <- seq(0, 1, length.out = bins + 1L)
  for (j in seq_len(d)) {
    edges <- unique(stats::quantile(xc[, j], probs, names = FALSE))
    cells[[j]] <- cut(xc[, j], edges, include.lowest = TRUE, labels = FALSE)
    midpoint[[j]] <- edges[-length(edges)] / 2 + edges[-1L] / 2
  }
  idx <- .irock_groups(c(cells, parts$discrete), n)
  representative <- vapply(idx, function(i) {
    center <- vapply(seq_len(d), function(j) midpoint[[j]][cells[[j]][i[1L]]],
                     numeric(1L))
    dx <- sweep(xc[i, , drop = FALSE], 2L, center, "-")
    # A common scale keeps the squared distance in range without changing
    # how coordinates are weighted.
    largest <- max(abs(dx))
    if (!is.finite(largest)) {
      stop("Bin distances overflowed; rescale the predictors.", call. = FALSE)
    }
    if (largest == 0) return(i[1L])
    i[which.min(rowSums((dx / largest)^2))]
  }, integer(1L))
  list(index = idx, representative = representative, continuous = xc,
       empirical = FALSE)
}

.irock_knn <- function(x, neighbors, centers, medoids, representative = NULL) {
  n <- nrow(x)
  if (is.null(neighbors)) neighbors <- min(n, 500L)
  if (is.null(centers)) centers <- min(n, 5L)
  if (neighbors > n || centers > n) {
    stop("'neighbors' and 'centers' cannot exceed the number of retained rows.",
         call. = FALSE)
  }
  xc <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  if (!ncol(xc)) stop("Nearest-neighbor fitting needs at least one predictor.",
                      call. = FALSE)
  if (as.double(neighbors) * centers > 5e6) {
    stop("Too many neighborhood memberships; reduce 'neighbors' or 'centers'.",
         call. = FALSE)
  }
  if (!is.null(representative)) {
    id <- representative
  } else if (centers == n) {
    id <- seq_len(n)
  } else if (medoids) {
    if (as.double(n) * (n - 1) / 2 > 1e7) {
      stop("The medoid distance matrix would be too large; use medoids = FALSE.",
           call. = FALSE)
    }
    if (!requireNamespace("fastkmedoids", quietly = TRUE)) {
      stop("Install 'fastkmedoids' to use medoids = TRUE.", call. = FALSE)
    }
    id <- 1L + attr(fastkmedoids::fastclarans(stats::dist(xc), n, centers),
                    "medoids")
  } else {
    id <- sample.int(n, centers)
  }
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Install 'FNN' to use binning = 'knn'.", call. = FALSE)
  }
  ii <- FNN::knnx.index(xc, xc[id, , drop = FALSE], k = neighbors)
  idx <- lapply(seq_along(id), function(j) {
    k <- ii[j, ]
    if (!id[j] %in% k) k[length(k)] <- id[j]
    k
  })
  list(index = idx, representative = id, continuous = xc, empirical = FALSE)
}
