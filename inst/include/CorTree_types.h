#ifndef CORTREE_TYPES_H
#define CORTREE_TYPES_H

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace cortree {

// Original determinant thresholds, evaluated on the log scale to avoid overflow.
inline bool update_precision_allowed(const arma::mat& precision) {
  double log_det = 0.0, sign = 0.0;
  arma::log_det(log_det, sign, precision);
  return !std::isfinite(log_det) || sign <= 0.0 || log_det < std::log(1e200);
}

inline bool regularize_precision(arma::mat& precision, arma::mat& covariance,
                                 double epsilon) {
  double log_det = 0.0, sign = 0.0;
  arma::log_det(log_det, sign, precision);
  if (sign <= 0.0 || !std::isfinite(log_det) || log_det <= std::log(1e150)) return false;
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, precision) || arma::any(eigenvalues <= 0))
    Rcpp::stop("Covariance regularization requires a positive-definite precision matrix.");
  arma::vec covariance_eigenvalues = 1.0 / eigenvalues + epsilon;
  precision = eigenvectors * arma::diagmat(1.0 / covariance_eigenvalues) * eigenvectors.t();
  covariance = eigenvectors * arma::diagmat(covariance_eigenvalues) * eigenvectors.t();
  return true;
}

inline double softplus(double x) {
  return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x)));
}

inline double binomial_log_kernel(double n, double y, double phi) {
  // Stable even when the logistic probability rounds to zero or one.
  return -y * softplus(-phi) - (n - y) * softplus(phi);
}

inline arma::vec independent_variance(const arma::mat& residual,
                                      const arma::vec& split_layer, double shape) {
  arma::vec rate = 1.0 / split_layer + arma::sum(arma::square(residual), 0).t() / 2.0;
  arma::vec out(rate.n_elem);
  for (arma::uword j = 0; j < rate.n_elem; ++j) {
    out(j) = 1.0 / R::rgamma(shape + residual.n_rows / 2.0, 1.0 / rate(j));
  }
  return out;
}

inline arma::vec horseshoe_nu(const arma::vec& lambda_sq) {
  arma::vec out(lambda_sq.n_elem);
  for (arma::uword j = 0; j < out.n_elem; ++j) {
    out(j) = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / lambda_sq(j)));
  }
  return out;
}

inline arma::vec stick_weights(const arma::uvec& counts, double alpha) {
  arma::vec pi(counts.n_elem, arma::fill::zeros);
  double remaining = 1.0;
  double tail = static_cast<double>(arma::accu(counts));
  for (arma::uword k = 0; k + 1 < counts.n_elem; ++k) {
    tail -= counts(k);
    double v = R::rbeta(1.0 + counts(k), alpha + tail);
    pi(k) = remaining * v;
    remaining *= 1.0 - v;
  }
  pi(counts.n_elem - 1) = remaining;
  return pi;
}

inline void validate_counts(const arma::mat& X) {
  if (X.n_rows == 0 || X.n_cols == 0 || !X.is_finite() ||
      arma::any(arma::vectorise(X) < 0) ||
      arma::any(arma::vectorise(X) != arma::floor(arma::vectorise(X)))) {
    Rcpp::stop("Counts must be a nonempty matrix of finite nonnegative integers.");
  }
  if (arma::any(arma::sum(X, 1) > std::numeric_limits<int>::max())) {
    Rcpp::stop("Row totals exceed the supported Polya-Gamma count range.");
  }
}

inline void validate_sampler(const arma::mat& X, int K, int total, int burnin,
                             int warm, int interval, double shape, double mu_var,
                             const arma::uvec& Z) {
  validate_counts(X);
  if (K < 1 || total < 1 || burnin < 0 || burnin >= total ||
      warm < 0 || warm > burnin || interval < 1) {
    Rcpp::stop("Require n_clus >= 1, 0 <= warm_start <= burnin < total_iter, and cov_interval >= 1.");
  }
  if (!std::isfinite(shape) || shape <= 0 || !std::isfinite(mu_var) || mu_var <= 0) {
    Rcpp::stop("Variance hyperparameters must be finite and positive.");
  }
  if (Z.n_elem != 1 && Z.n_elem != X.n_rows) {
    Rcpp::stop("init_Z must have one entry or nrow(X) entries.");
  }
  if (Z.n_elem == 0 || arma::any(Z >= static_cast<arma::uword>(K))) {
    Rcpp::stop("init_Z labels must be in 0, ..., n_clus - 1.");
  }
}

inline void validate_depth(int depth, int cutoff, bool all_ind) {
  if (depth < 1 || depth > 20 || (!all_ind && (cutoff < 0 || cutoff >= depth))) {
    Rcpp::stop("Require 1 <= tree_depth <= 20 and 0 <= cutoff_layer < tree_depth for correlated fits.");
  }
}
} // namespace cortree
#endif
