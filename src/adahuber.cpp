# include <RcppArmadillo.h>
# include <algorithm>
# include <string>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

// [[Rcpp::export]]
int sgn(const double x) {
  return (x > 0) - (x < 0);
}

// [[Rcpp::export]]
double f1(const double x, const arma::vec& resSq, const int n, const double rhs) {
  return arma::mean(arma::min(resSq / x, arma::ones(n))) - rhs;
}

// [[Rcpp::export]]
double rootf1(const arma::vec& resSq, const int n, const double rhs, double low, double up, const double tol = 0.001, const int maxIte = 500) {
  int ite = 1;
  while (ite <= maxIte && up - low > tol) {
    double mid = 0.5 * (up + low);
    double val = f1(mid, resSq, n, rhs);
    if (val < 0) {
      up = mid;
    } else {
      low = mid;
    }
    ite++;
  }
  return 0.5 * (low + up);
}

// [[Rcpp::export]]
double huberDer(const arma::vec& res, const double tau, const int n) {
  double rst = 0.0;
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    rst -= std::abs(cur) <= tau ? cur : tau * sgn(cur);
  }
  return rst / n;
}

// [[Rcpp::export]]
double huberMean(arma::vec X, const int n, const double tol = 0.0001, const int iteMax = 500) {
  double rhs = std::log(n) / n;
  double mx = arma::mean(X);
  X -= mx;
  double tau = arma::stddev(X) * std::sqrt((long double)n / std::log(n));
  double derOld = huberDer(X, tau, n);
  double mu = -derOld, muDiff = -derOld;
  arma::vec res = X - mu;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  double derNew = huberDer(res, tau, n);
  double derDiff = derNew - derOld;
  int ite = 1;
  while (std::abs(derNew) > tol && ite <= iteMax) {
    double alpha = 1.0;
    double cross = muDiff * derDiff;
    if (cross > 0) {
      double a1 = cross / derDiff * derDiff;
      double a2 = muDiff * muDiff / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    derOld = derNew;
    muDiff = -alpha * derNew;
    mu += muDiff;
    res = X - mu;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    derNew = huberDer(res, tau, n);
    derDiff = derNew - derOld;
    ite++;
  }
  return mu + mx;
}

// [[Rcpp::export]]
double mad(const arma::vec& x) {
  return 1.482602 * arma::median(arma::abs(x - arma::median(x)));
}

// [[Rcpp::export]]
arma::mat standardize(arma::mat X, const arma::rowvec& mx, const arma::vec& sx1, const int p) {
  for (int i = 0; i < p; i++) {
    X.col(i) = (X.col(i) - mx(i)) * sx1(i);
  }
  return X;
}

// [[Rcpp::export]]
void updateHuber(const arma::mat& Z, const arma::vec& res, arma::vec& der, arma::vec& grad, const int n, const double tau, const double n1) {
  for (int i = 0; i < n; i++) {
    double cur = res(i);
    if (std::abs(cur) <= tau) {
      der(i) = -cur;
    } else {
      der(i) = -tau * sgn(cur);
    }
  }
  grad = n1 * Z.t() * der;
}

// [[Rcpp::export]]
Rcpp::List adaHuberReg(const arma::mat& X, arma::vec Y, const double epsilon = 0.0001, const int iteMax = 500) {
  const int n = X.n_rows;
  const int p = X.n_cols;
  const double n1 = 1.0 / n;
  double rhs = n1 * (p + std::log(n * p));
  arma::rowvec mx = arma::mean(X, 0);
  arma::vec sx1 = 1.0 / arma::stddev(X, 0, 0).t();
  double my = arma::mean(Y);
  arma::mat Z = arma::join_rows(arma::ones(n), standardize(X, mx, sx1, p));
  Y -= my;
  double tau = 1.345 * mad(Y);
  arma::vec der(n);
  arma::vec gradOld(p + 1), gradNew(p + 1);
  updateHuber(Z, Y, der, gradOld, n, tau, n1);
  arma::vec beta = -gradOld, betaDiff = -gradOld;
  arma::vec res = Y - Z * beta;
  arma::vec resSq = arma::square(res);
  tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
  updateHuber(Z, res, der, gradNew, n, tau, n1);
  arma::vec gradDiff = gradNew - gradOld;
  int ite = 1;
  while (arma::norm(gradNew, "inf") > epsilon && ite <= iteMax) {
    double alpha = 1.0;
    double cross = arma::as_scalar(betaDiff.t() * gradDiff);
    if (cross > 0) {
      double a1 = cross / arma::as_scalar(gradDiff.t() * gradDiff);
      double a2 = arma::as_scalar(betaDiff.t() * betaDiff) / cross;
      alpha = std::min(std::min(a1, a2), 100.0);
    }
    gradOld = gradNew;
    betaDiff = -alpha * gradNew;
    beta += betaDiff;
    res -= Z * betaDiff;
    resSq = arma::square(res);
    tau = std::sqrt((long double)rootf1(resSq, n, rhs, arma::min(resSq), arma::accu(resSq)));
    updateHuber(Z, res, der, gradNew, n, tau, n1);
    gradDiff = gradNew - gradOld;
    ite++;
  }
  beta.rows(1, p) %= sx1;
  beta(0) = huberMean(Y + my - X * beta.rows(1, p), n, epsilon, iteMax);
  return Rcpp::List::create(Rcpp::Named("coef") = beta, Rcpp::Named("tau") = tau, Rcpp::Named("iteration") = ite);
}
