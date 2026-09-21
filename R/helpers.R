# Internal helper functions adapted from the esreg package.
.lesreg_conditional_truncated_variance <- function(y, x, approach) {
  if (sum(y <= 0) <= 2) {
    stop("Not enough negative quantile residuals!")
  }

  if (approach == "ind") {
    cv <- rep(stats::var(y[y <= 0]), length(y))
  } else {
    cv <- tryCatch({
      # Get conditional mean and sigma
      mu_sigma <- .lesreg_conditional_mean_sigma(y, x)
      mu <- mu_sigma$mu
      sigma <- mu_sigma$sigma

      # Truncated conditional variance
      if (approach == "scl_N") {
        beta <- -mu / sigma
        beta[beta < -30] <- -30
        cv <- sigma^2 * (1 - beta * stats::dnorm(beta)/stats::pnorm(beta) -
                           (stats::dnorm(beta)/stats::pnorm(beta))^2)
      } else if (approach == "scl_sp") {
        # Kernel density estimate of the standardized residuals
        df <- stats::density((y - mu) / sigma, bw = "SJ")
        lower_int_boundary <- min(df$x)
        pdf <- stats::approxfun(df, yleft = 0, yright = 0)

        # Truncation points
        beta <- -mu / sigma

        # Integrate the truncated pdf by the trapezoidal rule
        div <- 1000
        h <- (max(beta) - lower_int_boundary) / (div - 1)
        b_approx <- seq(lower_int_boundary, max(beta) + h, h)
        midpoint <- b_approx[-div] + h/2
        y0 <- pdf(b_approx)
        y1 <- b_approx * y0
        y2 <- b_approx^2 * y0
        cb <- cumsum(y0[-1] + y0[-div]) / 2 * h
        m1 <- cumsum(y1[-1] + y1[-div]) / 2 * h
        m2 <- cumsum(y2[-1] + y2[-div]) / 2 * h
        cv_approx <- m2 / cb - (m1 / cb)^2
        cv_approx[cv_approx < 0] <- NA

        # Approximate the conditional truncated variance
        cv <- sigma^2 * stats::approx(x = midpoint, y = cv_approx, xout = beta)$y
      }
      if (any(is.na(cv)) | any(!is.finite(cv)) | any(cv < 0)) stop() else cv
    }, error = function(e) {
      warning(paste0("Can not fit the ", approach, " estimator, switching to the ind approach!"))
      rep(stats::var(y[y <= 0]), length(y))
    })
  }

  cv
}

.lesreg_conditional_mean_sigma <- function(y, x) {
  # Starting values and ensure positive fitted standard deviations
  fit1 <- stats::lm(y ~ x - 1)
  fit2 <- stats::lm(abs(fit1$residuals) ~ x - 1)
  fit2$coefficients[1] <- fit2$coefficients[1] - min(0.001, min(fit2$fitted.values))
  b0 <- c(fit1$coefficients, fit2$coefficients)

  # Estimate the model under normality
  ll <- function(par, y, x) {
    k <- ncol(x)
    mu <- as.numeric(x %*% par[1:k])
    sigma <- as.numeric(x %*% par[(k+1):(2*k)])
    ifelse(all(sigma > 0), -sum(stats::dnorm(x=y, mean=mu, sd=sigma, log=TRUE)), NA)
  }
  fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="BFGS"), silent=TRUE)
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    fit <- try(stats::optim(b0, function(b) ll(par=b, y=y, x=x), method="Nelder-Mead",
                            control=list(maxit=10000)), silent=TRUE)
  }
  if(inherits(fit, "try-error") || (fit$convergence != 0)) {
    warning("Conditional variance estimation: cannot optimize model, returning starting values instead")
    b <- b0
  } else {
    b <- fit$par
  }

  # Estimated means and standard deviations
  k <- ncol(x)
  mu <- as.numeric(x %*% b[1:k])
  sigma <- as.numeric(x %*% b[(k+1):(2*k)])

  list(mu = mu, sigma = sigma)
}
