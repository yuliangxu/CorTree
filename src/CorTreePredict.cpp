#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

inline double logsumexp_vec(const arma::vec& x) {
  arma::uvec finite_idx = arma::find_finite(x);
  if (finite_idx.n_elem == 0) {
    return -arma::datum::inf;
  }
  double xmax = x(finite_idx).max();
  return xmax + std::log(arma::sum(arma::exp(x(finite_idx) - xmax)));
}

inline double logmeanexp_vec(const arma::vec& x) {
  double lse = logsumexp_vec(x);
  if (!arma::is_finite(lse)) {
    return -arma::datum::inf;
  }
  return lse - std::log(static_cast<double>(x.n_elem));
}

struct TreeSuffStats {
  arma::mat count;
  arma::mat kappa;
};

TreeSuffStats build_tree_suff_stats(const arma::mat& X, int tree_depth) {
  TreeSuffStats out;
  int total_nodes = static_cast<int>(std::pow(2.0, tree_depth + 1.0)) - 1;
  int total_parents = static_cast<int>(std::pow(2.0, tree_depth)) - 1;

  out.count = arma::zeros<arma::mat>(X.n_rows, static_cast<arma::uword>(total_nodes));
  out.kappa = arma::zeros<arma::mat>(X.n_rows, static_cast<arma::uword>(total_parents));

  arma::uvec interval_left(static_cast<arma::uword>(total_nodes), arma::fill::zeros);
  arma::uvec interval_right(static_cast<arma::uword>(total_nodes), arma::fill::zeros);

  int k = static_cast<int>(X.n_cols);
  out.count.col(0) = arma::sum(X, 1);
  interval_left(0) = 0;
  interval_right(0) = static_cast<arma::uword>(k);

  for (int i = 0; i < total_parents; i++) {
    arma::uword left_node = static_cast<arma::uword>(2 * i + 1);
    arma::uword right_node = static_cast<arma::uword>(2 * i + 2);

    interval_left(left_node) = interval_left(static_cast<arma::uword>(i));
    interval_right(left_node) = (interval_left(static_cast<arma::uword>(i)) + interval_right(static_cast<arma::uword>(i))) / 2;

    interval_left(right_node) = interval_right(left_node);
    interval_right(right_node) = interval_right(static_cast<arma::uword>(i));

    if (interval_left(static_cast<arma::uword>(i)) < interval_right(left_node)) {
      arma::uvec left_idx = arma::regspace<arma::uvec>(interval_left(static_cast<arma::uword>(i)), interval_right(left_node) - 1);
      out.count.col(left_node) = arma::sum(X.cols(left_idx), 1);
    }

    if (interval_left(right_node) < interval_right(static_cast<arma::uword>(i))) {
      arma::uvec right_idx = arma::regspace<arma::uvec>(interval_left(right_node), interval_right(static_cast<arma::uword>(i)) - 1);
      out.count.col(right_node) = arma::sum(X.cols(right_idx), 1);
    }

    out.kappa.col(static_cast<arma::uword>(i)) =
      out.count.col(left_node) - out.count.col(static_cast<arma::uword>(i)) / 2.0;
  }

  return out;
}

inline double log_binom_term(double n, double y) {
  return std::lgamma(n + 1.0) - std::lgamma(y + 1.0) - std::lgamma(n - y + 1.0);
}

double log_tree_lik_from_phi(const arma::rowvec& count_parent,
                             const arma::rowvec& kappa,
                             const arma::vec& phi) {
  double out = 0.0;
  for (arma::uword j = 0; j < phi.n_elem; ++j) {
    double n = count_parent(j);
    double y = kappa(j) + 0.5 * n;
    double phi_j = phi(j);
    double logp = -std::log1p(std::exp(-phi_j));
    double log1mp = -std::log1p(std::exp(phi_j));
    out += log_binom_term(n, y) + y * logp + (n - y) * log1mp;
  }
  return out;
}

arma::vec sample_phi_cluster(const arma::vec& mu_k,
                             const arma::mat& Sigma_inv_k,
                             const arma::vec& sigma2_k,
                             bool all_ind,
                             arma::uword L) {
  arma::vec phi(mu_k.n_elem, arma::fill::zeros);

  if (!all_ind && L > 0) {
    arma::mat chol_L = arma::chol(Sigma_inv_k, "lower");
    arma::vec z = arma::randn<arma::vec>(L);
    arma::vec noise = arma::solve(arma::trimatu(chol_L.t()), z, arma::solve_opts::fast);
    phi.subvec(0, L - 1) = mu_k.subvec(0, L - 1) + noise;
  }

  arma::uword ind_start = all_ind ? 0 : L;
  if (ind_start < mu_k.n_elem) {
    arma::uword ind_len = mu_k.n_elem - ind_start;
    arma::vec z_ind = arma::randn<arma::vec>(ind_len);
    phi.subvec(ind_start, mu_k.n_elem - 1) =
      mu_k.subvec(ind_start, mu_k.n_elem - 1) + arma::sqrt(sigma2_k) % z_ind;
  }

  phi.clamp(-7.0, 7.0);
  return phi;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List CorTree_heldout_logpred(const arma::mat& X_test,
                                   const Rcpp::List& mcmc,
                                   int tree_depth,
                                   int cutoff_layer,
                                   int burnin = 0,
                                   bool all_ind = false,
                                   int n_phi_mc = 1) {
  if (n_phi_mc <= 0) {
    Rcpp::stop("`n_phi_mc` must be positive.");
  }

  arma::cube mu = Rcpp::as<arma::cube>(mcmc["mu"]);
  arma::mat pi = Rcpp::as<arma::mat>(mcmc["pi"]);
  arma::cube sigma2_vec = Rcpp::as<arma::cube>(mcmc["sigma2_vec"]);

  arma::uword total_parents = mu.n_rows;
  arma::uword n_clus = mu.n_cols;
  arma::uword n_draw = mu.n_slices;
  arma::uword n_test = X_test.n_rows;
  arma::uword L = all_ind ? 0 : static_cast<arma::uword>(std::pow(2.0, cutoff_layer + 1.0) - 1.0);

  if (sigma2_vec.n_slices != n_draw) {
    Rcpp::stop("Mismatch: `sigma2_vec` and `mu` have different draw counts.");
  }
  if (pi.n_rows != n_clus) {
    Rcpp::stop("Mismatch: `pi` row count must equal n_clus.");
  }

  arma::uword pi_offset = 0;
  if (pi.n_cols == n_draw) {
    pi_offset = 0;
  } else if (pi.n_cols >= static_cast<arma::uword>(burnin) + n_draw) {
    pi_offset = static_cast<arma::uword>(burnin);
  } else {
    Rcpp::stop("`pi` does not have enough columns for requested burn-in offset.");
  }

  Rcpp::List Sigma_inv_sample;
  if (!all_ind) {
    Sigma_inv_sample = Rcpp::as<Rcpp::List>(mcmc["Sigma_inv"]);
    if (Sigma_inv_sample.size() < static_cast<int>(n_draw)) {
      Rcpp::stop("Mismatch: `Sigma_inv` list length is smaller than draw count.");
    }
  }

  TreeSuffStats test_stats = build_tree_suff_stats(X_test, tree_depth);
  arma::mat count_parent = test_stats.count.cols(0, total_parents - 1);
  arma::mat kappa = test_stats.kappa;

  arma::vec pointwise(n_test, arma::fill::zeros);

  for (arma::uword i = 0; i < n_test; ++i) {
    arma::vec log_draw(n_draw, arma::fill::zeros);
    arma::rowvec count_i = count_parent.row(i);
    arma::rowvec kappa_i = kappa.row(i);

    for (arma::uword s = 0; s < n_draw; ++s) {
      arma::mat mu_s = mu.slice(s);
      arma::mat sigma2_s = sigma2_vec.slice(s);
      arma::vec pi_s = pi.col(pi_offset + s);
      arma::vec log_mix(n_clus, arma::fill::zeros);

      arma::cube Sigma_inv_s;
      if (!all_ind) {
        Sigma_inv_s = Rcpp::as<arma::cube>(Sigma_inv_sample[static_cast<int>(s)]);
      }

      for (arma::uword k = 0; k < n_clus; ++k) {
        if (!arma::is_finite(pi_s(k)) || pi_s(k) <= 0.0) {
          log_mix(k) = -arma::datum::inf;
          continue;
        }

        arma::vec mu_k = mu_s.col(k);
        arma::vec sigma2_k;
        arma::mat Sigma_inv_k;

        if (all_ind) {
          sigma2_k = sigma2_s.col(k);
        } else {
          Sigma_inv_k = Sigma_inv_s.slice(k);
          sigma2_k = sigma2_s.col(k);
        }

        arma::vec log_phi_mc(static_cast<arma::uword>(n_phi_mc), arma::fill::zeros);
        for (int m = 0; m < n_phi_mc; ++m) {
          arma::vec phi_ik = sample_phi_cluster(mu_k, Sigma_inv_k, sigma2_k, all_ind, L);
          log_phi_mc(static_cast<arma::uword>(m)) = log_tree_lik_from_phi(count_i, kappa_i, phi_ik);
        }

        log_mix(k) = std::log(pi_s(k)) + logmeanexp_vec(log_phi_mc);
      }

      log_draw(s) = logsumexp_vec(log_mix);
    }

    pointwise(i) = logmeanexp_vec(log_draw);
  }

  return Rcpp::List::create(
    Rcpp::Named("hlpd") = arma::accu(pointwise),
    Rcpp::Named("pointwise_hlpd") = pointwise
  );
}

// [[Rcpp::export]]
Rcpp::List CorTree_heldout_membership(const arma::mat& X_test,
                                      const Rcpp::List& mcmc,
                                      int tree_depth,
                                      int cutoff_layer,
                                      int burnin = 0,
                                      bool all_ind = false,
                                      int n_phi_mc = 1) {
  if (n_phi_mc <= 0) {
    Rcpp::stop("`n_phi_mc` must be positive.");
  }

  arma::cube mu = Rcpp::as<arma::cube>(mcmc["mu"]);
  arma::mat pi = Rcpp::as<arma::mat>(mcmc["pi"]);
  arma::cube sigma2_vec = Rcpp::as<arma::cube>(mcmc["sigma2_vec"]);

  arma::uword total_parents = mu.n_rows;
  arma::uword n_clus = mu.n_cols;
  arma::uword n_draw = mu.n_slices;
  arma::uword n_test = X_test.n_rows;
  arma::uword L = all_ind ? 0 : static_cast<arma::uword>(std::pow(2.0, cutoff_layer + 1.0) - 1.0);

  if (sigma2_vec.n_slices != n_draw) {
    Rcpp::stop("Mismatch: `sigma2_vec` and `mu` have different draw counts.");
  }
  if (pi.n_rows != n_clus) {
    Rcpp::stop("Mismatch: `pi` row count must equal n_clus.");
  }

  arma::uword pi_offset = 0;
  if (pi.n_cols == n_draw) {
    pi_offset = 0;
  } else if (pi.n_cols >= static_cast<arma::uword>(burnin) + n_draw) {
    pi_offset = static_cast<arma::uword>(burnin);
  } else {
    Rcpp::stop("`pi` does not have enough columns for requested burn-in offset.");
  }

  Rcpp::List Sigma_inv_sample;
  if (!all_ind) {
    Sigma_inv_sample = Rcpp::as<Rcpp::List>(mcmc["Sigma_inv"]);
    if (Sigma_inv_sample.size() < static_cast<int>(n_draw)) {
      Rcpp::stop("Mismatch: `Sigma_inv` list length is smaller than draw count.");
    }
  }

  TreeSuffStats test_stats = build_tree_suff_stats(X_test, tree_depth);
  arma::mat count_parent = test_stats.count.cols(0, total_parents - 1);
  arma::mat kappa = test_stats.kappa;

  arma::vec pointwise(n_test, arma::fill::zeros);
  arma::mat post_prob(n_test, n_clus, arma::fill::zeros);
  arma::uvec z_hat(n_test, arma::fill::zeros);

  for (arma::uword i = 0; i < n_test; ++i) {
    arma::vec log_draw(n_draw, arma::fill::zeros);
    arma::mat draw_post(n_draw, n_clus, arma::fill::zeros);
    arma::rowvec count_i = count_parent.row(i);
    arma::rowvec kappa_i = kappa.row(i);

    for (arma::uword s = 0; s < n_draw; ++s) {
      arma::mat mu_s = mu.slice(s);
      arma::mat sigma2_s = sigma2_vec.slice(s);
      arma::vec pi_s = pi.col(pi_offset + s);
      arma::vec log_mix(n_clus, arma::fill::zeros);

      arma::cube Sigma_inv_s;
      if (!all_ind) {
        Sigma_inv_s = Rcpp::as<arma::cube>(Sigma_inv_sample[static_cast<int>(s)]);
      }

      for (arma::uword k = 0; k < n_clus; ++k) {
        if (!arma::is_finite(pi_s(k)) || pi_s(k) <= 0.0) {
          log_mix(k) = -arma::datum::inf;
          continue;
        }

        arma::vec mu_k = mu_s.col(k);
        arma::vec sigma2_k;
        arma::mat Sigma_inv_k;

        if (all_ind) {
          sigma2_k = sigma2_s.col(k);
        } else {
          Sigma_inv_k = Sigma_inv_s.slice(k);
          sigma2_k = sigma2_s.col(k);
        }

        arma::vec log_phi_mc(static_cast<arma::uword>(n_phi_mc), arma::fill::zeros);
        for (int m = 0; m < n_phi_mc; ++m) {
          arma::vec phi_ik = sample_phi_cluster(mu_k, Sigma_inv_k, sigma2_k, all_ind, L);
          log_phi_mc(static_cast<arma::uword>(m)) = log_tree_lik_from_phi(count_i, kappa_i, phi_ik);
        }

        log_mix(k) = std::log(pi_s(k)) + logmeanexp_vec(log_phi_mc);
      }

      double log_norm = logsumexp_vec(log_mix);
      log_draw(s) = log_norm;
      if (arma::is_finite(log_norm)) {
        draw_post.row(s) = arma::trans(arma::exp(log_mix - log_norm));
      }
    }

    pointwise(i) = logmeanexp_vec(log_draw);
    post_prob.row(i) = arma::mean(draw_post, 0);
    double row_sum = arma::accu(post_prob.row(i));
    if (row_sum > 0.0 && arma::is_finite(row_sum)) {
      post_prob.row(i) /= row_sum;
    } else {
      post_prob.row(i).fill(1.0 / static_cast<double>(n_clus));
    }
    z_hat(i) = post_prob.row(i).index_max();
  }

  return Rcpp::List::create(
    Rcpp::Named("hlpd") = arma::accu(pointwise),
    Rcpp::Named("pointwise_hlpd") = pointwise,
    Rcpp::Named("post_prob") = post_prob,
    Rcpp::Named("z_hat") = z_hat
  );
}
