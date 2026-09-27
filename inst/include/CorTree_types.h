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

inline void validate_ghs_diag_rate(double rate) {
  if (!std::isfinite(rate) || rate < 0.0)
    Rcpp::stop("ghs_diag_rate must be finite and nonnegative (0 selects legacy flat diagonals).");
}

inline void validate_ghs_diag_upper(double upper, double rate, int warm, bool all_ind) {
  if (std::isnan(upper) || upper <= 0.0)
    Rcpp::stop("ghs_diag_upper must be positive and finite, or Inf to disable the bound.");
  if (std::isfinite(upper) && rate > 0.0)
    Rcpp::stop("A finite ghs_diag_upper requires ghs_diag_rate = 0.");
  if (std::isfinite(upper) && !all_ind && warm > 0)
    Rcpp::stop("A finite ghs_diag_upper requires warm_start = 0 for correlated fits.");
}

inline double jmlr_ghs_tau_sq(double lambda) {
  double tau = 1.0 / lambda;
  return tau * tau;
}

inline void validate_ghs_jmlr_lambda(double lambda, double rate, double upper) {
  if (!std::isfinite(lambda) || lambda < 0.0)
    Rcpp::stop("ghs_jmlr_lambda must be finite and nonnegative (0 disables the JMLR prior).");
  if (lambda == 0.0) return;
  if (rate != 0.0 || std::isfinite(upper))
    Rcpp::stop("Positive ghs_jmlr_lambda requires ghs_diag_rate = 0 and ghs_diag_upper = Inf.");
  double tau_sq = jmlr_ghs_tau_sq(lambda);
  if (!(lambda / 2.0 > 0.0) || !(tau_sq > 0.0) || !std::isfinite(tau_sq) ||
      !std::isfinite(1.0 / tau_sq))
    Rcpp::stop("ghs_jmlr_lambda must give numerically representable diagonal and fixed global scales.");
}

inline double open_unit_uniform() {
  double u;
  do { u = R::runif(0.0, 1.0); } while (u <= 0.0 || u >= 1.0);
  return u;
}

inline void validate_ghs_det_df(double df, double upper, double rate = 0.0,
                              double fixed_scale = 0.0) {
  if (!std::isfinite(df) || df < 0.0)
    Rcpp::stop("ghs_det_df must be finite and nonnegative.");
  if (df > 0.0 && ((!std::isfinite(upper) && rate <= 0.0) || fixed_scale > 0.0))
    Rcpp::stop("Positive ghs_det_df requires a finite ghs_diag_upper or positive ghs_diag_rate, and no fixed global scale.");
}

inline void validate_ghs_scale_hierarchy(bool enabled, double shape,
                                        double rate_shape, double rate_rate,
                                        double diagonal_rate, double upper,
                                        double fixed_scale, int warm) {
  if (!std::isfinite(shape) || shape <= 0.0 ||
      !std::isfinite(rate_shape) || rate_shape <= 0.0 ||
      !std::isfinite(rate_rate) || rate_rate <= 0.0)
    Rcpp::stop("GHS scale hierarchy shape and rate hyperparameters must be finite and positive.");
  if (enabled && (diagonal_rate <= 0.0 || std::isfinite(upper) ||
                  fixed_scale != 0.0 || warm != 0))
    Rcpp::stop("ghs_scale_hierarchy requires positive ghs_diag_rate, infinite ghs_diag_upper, no fixed global scale, and warm_start = 0.");
}

inline double positive_gamma_draw(double shape, double rate) {
  double scale = 1.0 / rate;
  if (!std::isfinite(shape) || shape <= 0.0 || !std::isfinite(rate) || rate <= 0.0 ||
      !std::isfinite(scale) || scale <= 0.0)
    Rcpp::stop("GHS scale Gamma parameters must be finite, positive and representable.");
  double draw = R::rgamma(shape, scale);
  if (!std::isfinite(draw) || draw <= 0.0)
    Rcpp::stop("GHS scale Gamma draw must be finite and positive; no substitution was applied.");
  return draw;
}

// Omega = t Q. Scatter belongs to the effective Gaussian likelihood, without t.
inline double ghs_component_scale(const arma::mat& scatter,
                                  const arma::mat& template_precision, int n,
                                  double shape, double common_rate) {
  if (n < 0 || scatter.n_rows < 1 || scatter.n_cols != scatter.n_rows ||
      template_precision.n_rows != scatter.n_rows ||
      template_precision.n_cols != scatter.n_cols || !scatter.is_finite() ||
      !template_precision.is_finite() || arma::any(scatter.diag() < 0.0) ||
      !std::isfinite(shape) || shape <= 0.0 ||
      !std::isfinite(common_rate) || common_rate <= 0.0)
    Rcpp::stop("Invalid GHS component-scale conditional inputs.");
  double quadratic = arma::accu(scatter % template_precision.t());
  if (!std::isfinite(quadratic) || quadratic < 0.0)
    Rcpp::stop("GHS component-scale quadratic must be finite and nonnegative.");
  return positive_gamma_draw(shape + n * (scatter.n_rows / 2.0),
                             common_rate + quadratic / 2.0);
}

inline double ghs_common_scale_rate(const arma::vec& scales, double component_shape,
                                    double rate_shape, double rate_rate) {
  if (scales.n_elem < 1 || !scales.is_finite() || arma::any(scales <= 0.0) ||
      !std::isfinite(component_shape) || component_shape <= 0.0 ||
      !std::isfinite(rate_shape) || rate_shape <= 0.0 ||
      !std::isfinite(rate_rate) || rate_rate <= 0.0)
    Rcpp::stop("Invalid GHS shared scale-rate conditional inputs.");
  // Every instantiated component is included, even if its allocation is empty.
  return positive_gamma_draw(rate_shape + scales.n_elem * component_shape,
                             rate_rate + arma::accu(scales));
}

// Log integral of t^(shape-1) exp(-rate*t), up to constants independent of upper.
inline double bounded_gamma_log_mass(double shape, double rate, double upper) {
  if (!(upper > 0.0) || !std::isfinite(upper)) return -arma::datum::inf;
  if (rate == 0.0) return shape * std::log(upper);
  return R::pgamma(rate * upper, shape, 1.0, true, true);
}

inline double bounded_gamma_draw(double shape, double rate, double upper) {
  if (!(shape > 0.0) || !std::isfinite(shape) || rate < 0.0 ||
      !std::isfinite(rate) || !(upper > 0.0) || !std::isfinite(upper))
    Rcpp::stop("Invalid bounded GHS gamma parameters.");
  double log_u = std::log(open_unit_uniform());
  double value;
  if (rate == 0.0) {
    value = upper * std::exp(log_u / shape);
  } else {
    double log_mass = bounded_gamma_log_mass(shape, rate, upper);
    if (!std::isfinite(log_mass))
      Rcpp::stop("Bounded GHS gamma probability is not numerically representable.");
    value = R::qgamma(log_u + log_mass, shape, 1.0, true, true) / rate;
  }
  if (!std::isfinite(value) || value <= 0.0 || value >= upper)
    Rcpp::stop("Bounded GHS gamma draw is outside its open support; no clipping was applied.");
  return value;
}

struct BoundedGHSBlock {
  arma::vec beta;
  double gamma;
  arma::uword evaluations;
};

// Elliptical slice transition for beta after integrating out gamma, followed by
// its truncated-gamma conditional. This leaves the joint bounded block invariant.
inline BoundedGHSBlock bounded_ghs_block(const arma::vec& current,
                                        const arma::mat& inv_A,
                                        const arma::vec& inv_variance,
                                        const arma::vec& scatter_cross,
                                        double scatter_diagonal, double shape,
                                        double upper) {
  arma::mat precision = scatter_diagonal * inv_A + arma::diagmat(inv_variance);
  arma::mat factor;
  if (!precision.is_finite() || !arma::chol(factor, precision))
    Rcpp::stop("Bounded GHS beta conditional is not finite positive definite.");
  arma::vec rhs = arma::solve(arma::trimatl(factor.t()), scatter_cross, arma::solve_opts::fast);
  arma::vec mean = -arma::solve(arma::trimatu(factor), rhs, arma::solve_opts::fast);
  arma::vec direction = arma::solve(arma::trimatu(factor), arma::randn(current.n_elem), arma::solve_opts::fast);
  arma::vec centered = current - mean;
  double current_q = arma::dot(current, inv_A * current);
  double rate = scatter_diagonal / 2.0;
  double current_log_mass = bounded_gamma_log_mass(shape, rate, upper - current_q);
  if (!std::isfinite(current_q) || current_q < 0.0 || !std::isfinite(current_log_mass))
    Rcpp::stop("Bounded GHS current beta is outside its numerically representable support.");
  double log_slice = current_log_mass + std::log(open_unit_uniform());
  double angle = 2.0 * arma::datum::pi * open_unit_uniform();
  double lower_angle = angle - 2.0 * arma::datum::pi;
  double upper_angle = angle;
  for (arma::uword evaluations = 1; evaluations <= 10000; ++evaluations) {
    arma::vec beta = mean + centered * std::cos(angle) + direction * std::sin(angle);
    double q = arma::dot(beta, inv_A * beta);
    double available = upper - q;
    double log_mass = (std::isfinite(q) && q >= 0.0) ?
      bounded_gamma_log_mass(shape, rate, available) : -arma::datum::inf;
    if (std::isfinite(log_mass) && log_mass >= log_slice)
      return BoundedGHSBlock{beta, bounded_gamma_draw(shape, rate, available), evaluations};
    if (angle < 0.0) lower_angle = angle;
    else if (angle > 0.0) upper_angle = angle;
    else Rcpp::stop("Bounded GHS elliptical slice bracket stagnated.");
    angle = lower_angle + (upper_angle - lower_angle) * open_unit_uniform();
    if (evaluations % 64 == 0) Rcpp::checkUserInterrupt();
  }
  Rcpp::stop("Bounded GHS elliptical slice exceeded 10000 proposals; no state was substituted.");
}

inline void validate_ghs_precision_pair(const arma::mat& precision,
                                        const arma::mat& covariance) {
  arma::mat factor;
  if (!precision.is_finite() || !covariance.is_finite() ||
      !arma::chol(factor, precision) || !arma::chol(factor, covariance))
    Rcpp::stop("GHS precision/covariance is not finite positive definite; no regularization was applied.");
}

inline void validate_depth(int depth, int cutoff, bool all_ind) {
  if (depth < 1 || depth > 20 || (!all_ind && (cutoff < 0 || cutoff >= depth))) {
    Rcpp::stop("Require 1 <= tree_depth <= 20 and 0 <= cutoff_layer < tree_depth for correlated fits.");
  }
}
} // namespace cortree
#endif
