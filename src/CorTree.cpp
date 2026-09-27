#include "PolyaGamma.h" // to sample polya-gamma random variable
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "../inst/include/CorTree_types.h"



class CorTree{
private: 
  struct Data{
    arma::mat X; // n by (n_species) matrix, count table, sparse matrix
    int n_clus; // number of latent clusters 
    int n; // number of observations n
  } dat;
  
  struct FlatTree{
    int m; // tree depth. layer m=0 is the sample space Omega
    
    int cutoff_layer; // cutoff layer for correlated nodes, 0,...,m
    int L; // cutoff for the number of correlated nodes in a tree, any number in {2^{cutoff_layer+1}-1, i=0,1,...,m}
    arma::uvec idx_cor; // length of L
    arma::uvec idx_ind; // length of (total_parents-L)
    // for mu, fixed mu for the last layer
    arma::uvec idx_ind_mu; // length of (total_parents-L)
    arma::uvec idx_fixed_mu; // index of independent nodes in the last layer
    
    int total_nodes; // total number of nodes, = 2^{m+1}-1
    int total_parents; // total number of parents, = 2^m-1, also the total number of phi
    arma::mat count; // X_i(A_eps), n by total_nodes matrix, count table 
    arma::mat kappa; // n(A) - n(A_p)/2, n by total_parents matrix
    arma::uvec interval_left; // total_nodes vector, [l,r) for the interval of split on the node
    arma::uvec interval_right; // total_nodes vector, [l,r) for the interval of split on the node

    arma::vec ind_layer_idx;
    arma::vec cor_layer_idx;
  } tree; // fixed tree structure
  
  struct TreeLatentParameters{
    arma::mat mu; //  total_parents x n_clus, node probability
    arma::cube Sigma_inv; // L by L x n_clus precision matrix
    arma::mat sigma2_vec; // length of (total_parents-L) by n_clus matrix, diagonal elements of independent nodes. fixed
    arma::mat inv_a_vec;// for half-cauchy update
    arma::mat omega; // n by total_parents matrix, polya gamma latent variable
    
    // Dirichlet process  mixure parameters
    arma::uvec Z; // n by 1, membership indicator {1,2,...,n_clus}
    arma::uvec Z_cor_status; // n_clus by 1, whether to have Cov (1) or Independent (0)
    double alpha; // scalar, concentration parameter
    arma::vec pi; // n_clus by 1, mixture proportion
    arma::mat phi; // logit probability, n by total_parents matrix

    // trace log-like
    double loglike;
    arma::umat cluster_cor;
    arma::vec tau_sq_GHS;
  } paras;

  struct HyperParameters{
    double err_precision; // error precision for Sigma_inv
    double c_sigma2_vec; // decay constant for sigma2_vec
    double sigma_mu2; // prior variance for mu
  } hyper;

  
  struct GHS_params {
    arma::mat Omega;
    arma::mat Sigma;
    arma::mat Lambda_sq;
    arma::mat Nu;
    double tau_sq;
    double xi;
    
    
    // Constructor that takes 'L' (tree.L) and initializes each member accordingly
    GHS_params(int L) {
      Omega = arma::eye<arma::mat>(L, L);        // Identity matrix for Omega
      Sigma = arma::eye<arma::mat>(L, L);        // Identity matrix for Sigma
      Lambda_sq = arma::ones<arma::mat>(L, L);   // Ones matrix for Lambda_sq
      Nu = arma::ones<arma::mat>(L, L);          // Ones matrix for Nu
      tau_sq = 0.1;                              // Initialize tau_sq to 1.0
      xi = 1.0;                                  // Initialize xi to 1.0
    }
  };

  struct Gibbs_control{
    int total_iter;
    int burnin;
    int mcmc_sample;
  }gibbs_control;
  
  struct MCMCsample{
    arma::cube mu_sample;
    arma::cube phi_sample;
    arma::mat pi_sample;
    arma::cube sigma2_vec_sample;
    arma::umat Z_sample; 
    Rcpp::List Sigma_inv_sample;
    arma::vec loglik;
    arma::ucube cluster_cor;
    arma::mat tausq_GHS;
    arma::mat scale_t;
    arma::vec scale_b;
  } paras_sample;
  
public:
  
  arma::uvec covariance_regularizations;
  arma::uvec covariance_updates_skipped;
  arma::uvec empty_precision_updates;
  int iter=0;
  std::vector<GHS_params> ghs_list;
  std::vector<arma::uvec> ghs_ind_all;
  Rcpp::List test_output;
  bool all_ind = false;
  bool save_phi_trace = true;
  bool save_cluster_cor_trace = true;
  int warm_start = 0;
  int cov_interval = 1;
  double ghs_diag_rate = 0.0;
  double ghs_diag_upper = arma::datum::inf;
  double ghs_jmlr_lambda = 0.0;
  double ghs_det_df = 0.0;
  bool ghs_scale_hierarchy = false;
  double ghs_scale_shape = 3.0;
  double ghs_scale_rate_shape = 2.0;
  double ghs_scale_rate_rate = 1.0;
  arma::vec precision_scales;
  double precision_scale_rate = 2.0;
  arma::uword uniform_block_updates = 0;
  arma::uword uniform_ess_evaluations = 0;
  arma::uword uniform_max_ess_evaluations = 0;
  double z_det_gamma0 = 1.0;
  bool z_det_gamma_linear = false;
  int z_mode = 0; // 0: current robust update_Z, 1: v1-like update_Z
  int clustering_freq = 1;

  arma::uword find_parent(arma::uword node){
    return floor((node-1)/2);
  };
  
  arma::uvec find_parent(arma::uvec node){
    return arma::floor((node-1)/2);
  };
  
  arma::uword find_left_child(arma::uword node){
    return 2*node+1;
  };
  
  arma::uvec find_left_child(arma::uvec node){
    return 2*node+1;
  };
  
  arma::uword find_right_child(arma::uword node){
    return 2*node+2;
  };
  
  void load_data(arma::mat X, int n_clus){
    dat.X = X;
    dat.n_clus = n_clus;
    covariance_regularizations.zeros(n_clus);
    covariance_updates_skipped.zeros(n_clus);
    empty_precision_updates.zeros(n_clus);
    dat.n = X.n_rows;
  };

  void set_hyperparameter(double c_sigma2_vec, double sigma_mu2){
    hyper.c_sigma2_vec = c_sigma2_vec;
    hyper.sigma_mu2 = sigma_mu2;
    
  };

  void set_gibbs_control(int total_iter, int burnin){
    gibbs_control.total_iter = total_iter;
    gibbs_control.burnin = burnin;
    gibbs_control.mcmc_sample = total_iter - burnin;
  };

  void initialize_ghs_indices(int p){
    ghs_ind_all.clear();
    if(p <= 0){
      return;
    }
    ghs_ind_all.resize(static_cast<std::size_t>(p));
    for(int i = 0; i < p; ++i){
      arma::uvec ind = arma::regspace<arma::uvec>(0, p - 1);
      ind.shed_row(static_cast<arma::uword>(i));
      ghs_ind_all[static_cast<std::size_t>(i)] = std::move(ind);
    }
  }

  void initialize_mcmc_sample(){
    paras_sample.mu_sample = arma::zeros<arma::cube>(tree.total_parents, dat.n_clus, gibbs_control.mcmc_sample);
    if(save_phi_trace){
      paras_sample.phi_sample = arma::zeros<arma::cube>(dat.n, tree.total_parents, gibbs_control.mcmc_sample);
    }
    if(!all_ind){
      paras_sample.Sigma_inv_sample = Rcpp::List(gibbs_control.mcmc_sample);
    }
    // Save full pi trajectory by default (including burn-in iterations).
    paras_sample.pi_sample = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.total_iter);
    paras_sample.loglik = arma::zeros<arma::vec>(gibbs_control.total_iter);
    paras_sample.Z_sample = arma::zeros<arma::umat>(dat.n, gibbs_control.mcmc_sample);
    
    if(save_cluster_cor_trace){
      paras_sample.cluster_cor = arma::zeros<arma::ucube>(dat.n, dat.n, gibbs_control.mcmc_sample);
    }
    paras_sample.tausq_GHS = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.mcmc_sample);
    if (!all_ind && ghs_scale_hierarchy) {
      paras_sample.scale_t.zeros(dat.n_clus, gibbs_control.mcmc_sample);
      paras_sample.scale_b.zeros(gibbs_control.mcmc_sample);
    }
  };
  
  
  void vectorize_tree(arma::uword tree_depth, arma::uword cutoff_layer){
    // set tree structure
    tree.m = tree_depth;
    tree.total_nodes = pow(2, tree.m+1) - 1;
    tree.total_parents = pow(2, tree.m) - 1;
    
    tree.cutoff_layer = cutoff_layer;
    if(all_ind){
      tree.L = 0;
      tree.idx_ind = arma::regspace<arma::uvec>(0, tree.total_parents-1);
    }else{
      tree.L = pow(2, tree.cutoff_layer+1) - 1;
      tree.idx_cor = arma::regspace<arma::uvec>(0, tree.L-1);
      if (tree.L < tree.total_parents)
        tree.idx_ind = arma::regspace<arma::uvec>(tree.L, tree.total_parents-1);
      else tree.idx_ind.reset();
      initialize_ghs_indices(tree.L);
    }
    

    // ser precision for Sigma_inv
    hyper.err_precision = pow(10.0, -15.0);
    Rcout<<"err_precision="<<hyper.err_precision<<std::endl;
    
    // initialize the tree
    tree.count = arma::zeros<arma::mat>(dat.X.n_rows, tree.total_nodes);
    tree.kappa = arma::zeros<arma::mat>(dat.X.n_rows, tree.total_parents);
    // Old wrong code (fragile in this environment): this uvec init path has
    // caused runtime layout/init failures for some builds.
    // tree.interval_left = arma::zeros<arma::uvec>(tree.total_nodes);
    // tree.interval_right = arma::zeros<arma::uvec>(tree.total_nodes);
    tree.interval_left.set_size(static_cast<arma::uword>(tree.total_nodes));
    tree.interval_left.zeros();
    tree.interval_right.set_size(static_cast<arma::uword>(tree.total_nodes));
    tree.interval_right.zeros();
    
    // fill in the tree by the count matrix
    int k = dat.X.n_cols;
    tree.count.col(0) = arma::sum(dat.X, 1);
    tree.interval_left(0) = 0; // [left, right) left closed, right open
    tree.interval_right(0) = k;

    for(int i=0; i<tree.total_parents; i++){
      
      arma::uword left_node = find_left_child(i);
      tree.interval_left(left_node) = tree.interval_left(i);
      tree.interval_right(left_node) = (tree.interval_left(i) + tree.interval_right(i))/2;
      
      arma::uword right_node = find_right_child(i);
      tree.interval_left(right_node) = tree.interval_right(left_node);
      tree.interval_right(right_node) = tree.interval_right(i);
      
      // correct
      arma::uvec left_idx, right_idx;
      // Old wrong code: can underflow with unsigned arithmetic when right bound is 0.
      // if(tree.interval_left(i) <= tree.interval_right(left_node)-1){
      if(tree.interval_left(i) < tree.interval_right(left_node)){
        left_idx = arma::regspace<arma::uvec>(tree.interval_left(i), tree.interval_right(left_node)-1);
        tree.count.col(left_node) = arma::sum(dat.X.cols(left_idx), 1);
      }else{
        tree.count.col(left_node) = arma::zeros<arma::vec>(dat.X.n_rows);
      }
      
      
      // Old wrong code: can underflow with unsigned arithmetic when right bound is 0.
      // if(tree.interval_left(right_node) <= tree.interval_right(i)-1){
      if(tree.interval_left(right_node) < tree.interval_right(i)){
        right_idx = arma::regspace<arma::uvec>(tree.interval_left(right_node), tree.interval_right(i)-1);
        tree.count.col(right_node) = arma::sum(dat.X.cols(right_idx), 1);
      }else{
        tree.count.col(right_node) = arma::zeros<arma::vec>(dat.X.n_rows);
      }
      
      //wrong: did not consider overflow of nodes
      // arma::uvec left_idx = arma::regspace<arma::uvec>(tree.interval_left(i), tree.interval_right(left_node)-1);
      // arma::uvec right_idx = arma::regspace<arma::uvec>(tree.interval_left(right_node), tree.interval_right(i)-1);
      
      // tree.count.col(left_node) = arma::sum(dat.X.cols(left_idx), 1);
      // tree.count.col(right_node) = arma::sum(dat.X.cols(right_idx), 1);
      
      // get kappa
      tree.kappa.col(i) = tree.count.col(left_node) - tree.count.col(i)/2;

    }
    // Split probabilities belong to the child layer: root split is layer 1.
    arma::vec split_layer = 1.0 + arma::floor(arma::log2(
      1.0 + arma::regspace<arma::vec>(0, tree.total_parents - 1)));
    tree.ind_layer_idx = split_layer.elem(tree.idx_ind);
    tree.cor_layer_idx = split_layer.elem(tree.idx_cor);

  }; 

  void set_test_output(){
    test_output = Rcpp::List::create(Rcpp::Named("count") = tree.count,
                                     Rcpp::Named("sigma_vec") = paras.sigma2_vec,
                                     Rcpp::Named("kappa") = tree.kappa,
                                     Rcpp::Named("interval_left") = tree.interval_left,
                                     Rcpp::Named("interval_right") = tree.interval_right);
  }
  
  
  void initialize_parameters(arma::uvec init_Z){
    if (!all_ind && ghs_scale_hierarchy) {
      precision_scales.ones(dat.n_clus);
      precision_scale_rate = ghs_scale_shape > 1.0 ? ghs_scale_shape - 1.0 : ghs_scale_shape;
    }
    paras.mu = arma::zeros<arma::mat>(tree.total_parents, dat.n_clus);
    paras.Z_cor_status = arma::ones<arma::uvec>(dat.n_clus); // whether to have Cov or Independent
    // Assuming tree.L is the dimension of the identity matrix and tree.total_slices is the number of slices

    if(all_ind){
      paras.sigma2_vec = arma::ones<arma::mat>(tree.total_parents, dat.n_clus);
      paras.inv_a_vec = arma::ones<arma::mat>(tree.total_parents, dat.n_clus);
      paras_sample.sigma2_vec_sample = arma::zeros<arma::cube>(tree.total_parents, dat.n_clus, gibbs_control.mcmc_sample);
    }else{
      paras.Sigma_inv = arma::cube(tree.L, tree.L, dat.n_clus, arma::fill::zeros);
      for (size_t i = 0; i < dat.n_clus; i++) {
          paras.Sigma_inv.slice(i) = arma::eye<arma::mat>(tree.L, tree.L);
      }
      paras.sigma2_vec = arma::ones<arma::mat>(tree.total_parents - tree.L, dat.n_clus);
      paras.inv_a_vec = arma::ones<arma::mat>(tree.total_parents - tree.L, dat.n_clus);
      paras_sample.sigma2_vec_sample = arma::zeros<arma::cube>(tree.total_parents - tree.L, dat.n_clus, gibbs_control.mcmc_sample);
    }
    

    

    
    paras.omega = arma::zeros<arma::mat>(dat.X.n_rows, tree.total_parents);
    // paras.Z = arma::randi<arma::uvec>(dat.n, arma::distr_param(0, dat.n_clus-1)); // random initialization of membership
    // paras.Z = arma::randi<arma::uvec>(dat.X.n_rows, arma::distr_param(0, dat.n_clus - 1));
    paras.Z = init_Z;
    
    
    paras.alpha = 1.0; //
    paras.pi = arma::ones<arma::vec>(dat.n_clus)/dat.n_clus; //
    paras.phi = arma::zeros<arma::mat>(dat.n, tree.total_parents);


    for (int i = 0; i < dat.n_clus; i++) {
      ghs_list.push_back(GHS_params(tree.L));
      if (!all_ind && ghs_jmlr_lambda > 0.0) {
        // The JMLR common scale fixes tau=1/lambda; it is not sampled.
        ghs_list.back().tau_sq = cortree::jmlr_ghs_tau_sq(ghs_jmlr_lambda);
      }
      if (!all_ind && std::isfinite(ghs_diag_upper)) {
        double initial_precision = std::min(1.0, ghs_diag_upper / 2.0);
        if (initial_precision <= 0.0 || !std::isfinite(1.0 / initial_precision))
          Rcpp::stop("ghs_diag_upper is too small for finite precision/covariance initialization.");
        ghs_list.back().Omega *= initial_precision;
        ghs_list.back().Sigma /= initial_precision;
        paras.Sigma_inv.slice(i) = ghs_list.back().Omega;
      }
    }
    
    paras.tau_sq_GHS = arma::zeros<arma::vec>(dat.n_clus);
    if (!all_ind && ghs_jmlr_lambda > 0.0)
      paras.tau_sq_GHS.fill(cortree::jmlr_ghs_tau_sq(ghs_jmlr_lambda));

  };
  
  
  
  
  void update_omega(){
    
    // #pragma omp parallel for
    for(int i = 0; i < dat.X.n_rows; i++){
      arma::vec b_i = tree.count.row(i).cols(0, tree.total_parents-1).t();
      arma::vec c_i = paras.phi.row(i).t();
      arma::vec omega_i = arma_pgdraw(b_i, c_i);
      paras.omega.row(i) = omega_i.t();
    }
  };
  
  
  void update_phi(){

    // #pragma omp parallel for
    for(arma::uword i=0; i<dat.n; i++){
      arma::uword k = paras.Z(i);
      arma::vec omega_i = paras.omega.row(i).t();
      arma::vec kappa_i = tree.kappa.row(i).t();
      arma::vec phi_i = arma::zeros<arma::vec>(tree.total_parents);
      // ---------- correlated ---------- //
      // Update the diagonal of paras.Sigma_inv.slice(k) directly, without creating a full diagonal matrix
      if(!all_ind){
        arma::mat post_Sigma_inv = paras.Sigma_inv.slice(k);
        post_Sigma_inv.diag() += omega_i.elem(tree.idx_cor);  // In-place addition to the diagonal


        arma::mat L = arma::chol(post_Sigma_inv, "lower");
        // Compute the right-hand side vector 'b'
        arma::vec b = paras.Sigma_inv.slice(k) * paras.mu(tree.idx_cor, arma::uvec{k}) + kappa_i.elem(tree.idx_cor);
        // Solve L * y = b and then L.t() * post_mu = y
        arma::vec y = arma::solve(arma::trimatl(L), b, arma::solve_opts::fast);      // Solve L * y = b
        arma::vec post_mu = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);  // Solve L.t() * post_mu = y
        // Generate the random vector and solve in one step
        arma::vec z = arma::randn<arma::vec>(tree.L);
        arma::vec random_vec = arma::solve(arma::trimatu(L.t()), z, arma::solve_opts::fast);  // Solve L.t() * random_vec = z
        // Update phi_i in place
        phi_i(tree.idx_cor) = post_mu + random_vec;
      }
      

      

      // ---------- independent ---------- //
      arma::vec post_sigma2 = 1/(1/paras.sigma2_vec.col(k) + omega_i.elem(tree.idx_ind));
      arma::vec post_mu2 = post_sigma2 % (paras.mu(tree.idx_ind, arma::uvec{k})/paras.sigma2_vec.col(k) + 
        kappa_i.elem(tree.idx_ind));
      phi_i(tree.idx_ind) = post_mu2 +
        sqrt(post_sigma2) % arma::randn<arma::vec>(post_mu2.n_elem);

      // Draw from the Gaussian full conditional without clipping.

      paras.phi.row(i) = phi_i.t();
    }
    
  };
  
  void update_mu(){
    // #pragma omp parallel for
    for(arma::uword k=0; k<dat.n_clus; k++){
      arma::uvec idx_k = arma::find(paras.Z == k);
      {
        // correlated
        if(!all_ind){
          arma::mat post_Sigma_inv = idx_k.n_elem * paras.Sigma_inv.slice(k);
          post_Sigma_inv.diag() += 1/hyper.sigma_mu2;
          arma::mat L = arma::chol(post_Sigma_inv, "lower");
          arma::vec b = paras.Sigma_inv.slice(k)*arma::sum(paras.phi(idx_k,tree.idx_cor), 0).t();
          arma::vec y = arma::solve(arma::trimatl(L), b, arma::solve_opts::fast);
          arma::vec post_mu = arma::solve(arma::trimatu(L.t()), y, arma::solve_opts::fast);
          arma::vec z = arma::randn<arma::vec>(tree.idx_cor.n_elem);
          arma::vec random_vec = arma::solve(arma::trimatu(L.t()), z, arma::solve_opts::fast);
          paras.mu(tree.idx_cor, arma::uvec{k}) = post_mu + random_vec;
        }

        
        // independent
        arma::vec post_sigma2 = 1/(idx_k.n_elem/paras.sigma2_vec.col(k) + 1/hyper.sigma_mu2);
        arma::vec post_mu2 = post_sigma2 % (arma::sum(paras.phi(idx_k,tree.idx_ind), 0).t()/paras.sigma2_vec.col(k));
        paras.mu(tree.idx_ind, arma::uvec{k}) = post_mu2 +
          sqrt(post_sigma2) % arma::randn<arma::vec>(tree.idx_ind.n_elem);

      }
      
    }
  }; // mu is set to 0 for now

  void GHS_oneSample(const arma::mat& S, int n,
                   GHS_params& params) {
    int p = S.n_rows;
    if (ghs_jmlr_lambda > 0.0)
      params.tau_sq = cortree::jmlr_ghs_tau_sq(ghs_jmlr_lambda);
    if (n < 0 || p < 1 || S.n_cols != S.n_rows || !S.is_finite() ||
        arma::any(S.diag() < 0.0))
      Rcpp::stop("GHS update requires a finite scatter matrix with nonnegative diagonal.");
    bool bounded = std::isfinite(ghs_diag_upper);
    if (bounded && (arma::any(params.Omega.diag() <= 0.0) ||
                    arma::any(params.Omega.diag() >= ghs_diag_upper)))
      Rcpp::stop("Bounded GHS precision diagonal is outside (0, ghs_diag_upper).");
    arma::vec effective_diagonal = S.diag() + 2.0 * ghs_diag_rate;
    if (!effective_diagonal.is_finite() || (!bounded &&
        (arma::any(effective_diagonal <= 0.0) || (ghs_diag_rate == 0.0 && n < 1))))
      Rcpp::stop("GHS update requires positive diagonal rates; an empty cluster requires ghs_diag_rate > 0.");
    if (p == 1) {
      params.Omega(0,0) = bounded ?
        cortree::bounded_gamma_draw((n + ghs_det_df) / 2.0 + 1.0, S(0,0) / 2.0, ghs_diag_upper) :
        R::rgamma((n + ghs_det_df) / 2.0 + 1.0, 2.0 / effective_diagonal(0));
      if (!std::isfinite(params.Omega(0,0)) || params.Omega(0,0) <= 0.0)
        Rcpp::stop("GHS precision draw must be finite and positive.");
      params.Sigma(0,0) = 1.0 / params.Omega(0,0);
      if (ghs_diag_rate > 0.0 || bounded)
        cortree::validate_ghs_precision_pair(params.Omega, params.Sigma);
      return;
    }

    if(static_cast<int>(ghs_ind_all.size()) != p){
      initialize_ghs_indices(p);
    }
    
    // Sample Sigma and Omega = inv(Sigma)
    for (arma::uword i = 0; i < p; ++i) {
      const arma::uvec& ind = ghs_ind_all[static_cast<std::size_t>(i)];
      
      // Extract submatrices using proper matrix subsetting
      arma::mat Sigma_11 = params.Sigma(ind, ind);             // Sigma_11
      arma::vec sigma_12 = params.Sigma(ind, arma::uvec{i});    // Sigma_12
      double sigma_22 = params.Sigma(i,i);  // Sigma_22 as submatrix
      
      arma::vec s_21 = S(ind, arma::uvec{i});            // s_21
      // exp(-rho * tr(Omega)) shifts both the Schur gamma rate and beta precision.
      double s_22 = effective_diagonal(i);
      
      arma::vec lambda_sq_12 = params.Lambda_sq(ind, arma::uvec{i});    // Lambda_sq_12
      arma::vec nu_12 = params.Nu(ind, arma::uvec{i});                  // Nu_12
      
      double gamma;
      arma::mat inv_Omega_11;
      arma::vec beta;
      if (bounded) {
        arma::mat Omega_11 = params.Omega(ind, ind);
        if (!arma::inv_sympd(inv_Omega_11, Omega_11))
          Rcpp::stop("Bounded GHS precision subblock inversion failed.");
        arma::vec current_beta = params.Omega(ind, arma::uvec{i});
        cortree::BoundedGHSBlock draw = cortree::bounded_ghs_block(
          current_beta, inv_Omega_11, 1.0 / (lambda_sq_12 * params.tau_sq),
          s_21, S(i,i), (n + ghs_det_df) / 2.0 + 1.0, ghs_diag_upper);
        beta = draw.beta;
        gamma = draw.gamma;
        ++uniform_block_updates;
        uniform_ess_evaluations += draw.evaluations;
        uniform_max_ess_evaluations = std::max(uniform_max_ess_evaluations, draw.evaluations);
      } else {
        // Sample gamma and beta using arma::randg
        gamma = arma::randg<double>(arma::distr_param((n + ghs_det_df) / 2.0 + 1, 2.0 / s_22));
        if (!std::isfinite(gamma) || gamma <= 0.0)
          Rcpp::stop("GHS Schur-complement precision draw must be finite and positive.");
        inv_Omega_11 = Sigma_11 - sigma_12 * sigma_12.t() / sigma_22;  // Convert sigma_22 to scalar

        arma::mat inv_C = s_22 * inv_Omega_11 + arma::diagmat(1.0 / (lambda_sq_12 * params.tau_sq));
        arma::mat inv_C_chol = arma::chol(inv_C);
        // Reuse the Cholesky factor for the mean as well as the noise.
        // A generic solve can silently substitute a pseudoinverse for an
        // ill-conditioned SPD matrix, changing this Gaussian conditional.
        arma::vec rhs = arma::solve(arma::trimatl(inv_C_chol.t()), s_21, arma::solve_opts::fast);
        arma::vec mu_i = -arma::solve(arma::trimatu(inv_C_chol), rhs, arma::solve_opts::fast);
        beta = mu_i + arma::solve(arma::trimatu(inv_C_chol), arma::randn(p - 1), arma::solve_opts::fast);
      }
      arma::vec omega_12 = beta;
      double omega_22 = gamma + arma::as_scalar(beta.t() * inv_Omega_11 * beta);
      if (bounded && (!std::isfinite(omega_22) || omega_22 >= ghs_diag_upper))
        Rcpp::stop("Bounded GHS diagonal draw is outside its support; no clipping was applied.");
      
      // Update Lambda_sq and Nu using arma::randg
      arma::vec rate = arma::square(omega_12) / (2.0 * params.tau_sq) + 1.0 / nu_12;
      arma::mat inv_lambda_sq_12 = arma::randg<arma::vec>(rate.n_elem, arma::distr_param(1.0,1.0)); 
      inv_lambda_sq_12 /= rate;
      lambda_sq_12 = 1.0 / inv_lambda_sq_12;
      nu_12 = cortree::horseshoe_nu(lambda_sq_12);

      // Store omega_12 and omega_22 in Omega matrix
      params.Omega(arma::uvec{i}, ind) = omega_12.t();
      params.Omega(ind, arma::uvec{i}) = omega_12;
      params.Omega(i, i) = omega_22;
      
      
      
      // Update Sigma
      arma::vec temp = inv_Omega_11*beta;
      Sigma_11 = inv_Omega_11 + temp*temp.t()/gamma;
      sigma_12 = -temp/gamma; 
      sigma_22 = 1/gamma;
      params.Sigma(ind,ind) = Sigma_11; 
      params.Sigma(i,i) = sigma_22;
      params.Sigma(arma::uvec{i},ind) = sigma_12.t(); 
      params.Sigma(ind,arma::uvec{i}) = sigma_12;
      
      // Update Lambda_sq matrix
      params.Lambda_sq(arma::uvec{i},ind) = lambda_sq_12.t();
      params.Lambda_sq(ind, arma::uvec{i}) = lambda_sq_12;
      params.Nu(arma::uvec{i},ind) = nu_12.t(); 
      params.Nu(ind,arma::uvec{i}) = nu_12;
      
      // The JMLR prior has a fixed common scale, so no tau or xi draw.
      if (ghs_jmlr_lambda == 0.0) {
        arma::vec omega_vector = params.Omega(arma::trimatl_ind(size(params.Omega), -1));
        arma::vec lambda_sq_vector = params.Lambda_sq(arma::trimatl_ind(size(params.Lambda_sq), -1));
        double rate_tau_sq = 1.0 / params.xi + arma::sum(arma::square(omega_vector) / (2.0 * lambda_sq_vector));
        params.tau_sq = 1.0 / arma::randg<double>(arma::distr_param((p * (p - 1) / 2 + 1) / 2.0, 1.0 / rate_tau_sq));
        params.xi = 1.0 / arma::randg<double>(arma::distr_param(1.0, 1.0 / (1.0 + 1.0 / params.tau_sq)));
      }
    }

    
    if (ghs_diag_rate > 0.0 || bounded)
      cortree::validate_ghs_precision_pair(params.Omega, params.Sigma);
}

  
  void update_Sigma(){
    for (arma::uword k = 0; k < static_cast<arma::uword>(dat.n_clus); ++k) {
      arma::uvec idx_k = arma::find(paras.Z == k);
      if (!all_ind && (idx_k.n_elem > 0 || ghs_diag_rate > 0.0 || std::isfinite(ghs_diag_upper))) {
        arma::mat residual = paras.phi(idx_k, tree.idx_cor);
        residual.each_row() -= paras.mu(tree.idx_cor, arma::uvec{k}).t();
        // The diagonal warm start is burn-in initialization, not a GHS update.
        if (iter < warm_start && idx_k.n_elem > 0) {
          arma::vec variance = cortree::independent_variance(
            residual, tree.cor_layer_idx, hyper.c_sigma2_vec);
          paras.Sigma_inv.slice(k) = arma::diagmat(1.0 / variance);
          // Keep the block sampler's inverse pair consistent across warm-up.
          ghs_list[k].Omega = paras.Sigma_inv.slice(k);
          ghs_list[k].Sigma = arma::diagmat(variance);
        } else if (ghs_diag_rate > 0.0 || std::isfinite(ghs_diag_upper) ||
                   cortree::update_precision_allowed(paras.Sigma_inv.slice(k))) {
          arma::mat scatter = residual.t() * residual;
          if (ghs_scale_hierarchy) {
            // The GHS state remains Q and Q^{-1}; likelihoods use Omega = t Q.
            GHS_oneSample(precision_scales(k) * scatter,
                          static_cast<int>(idx_k.n_elem), ghs_list[k]);
            precision_scales(k) = cortree::ghs_component_scale(
              scatter, ghs_list[k].Omega, static_cast<int>(idx_k.n_elem),
              ghs_scale_shape, precision_scale_rate);
            paras.Sigma_inv.slice(k) = precision_scales(k) * ghs_list[k].Omega;
            arma::mat effective_covariance = ghs_list[k].Sigma / precision_scales(k);
            cortree::validate_ghs_precision_pair(paras.Sigma_inv.slice(k), effective_covariance);
          } else {
            GHS_oneSample(scatter, static_cast<int>(idx_k.n_elem), ghs_list[k]);
            paras.Sigma_inv.slice(k) = ghs_list[k].Omega;
          }
          paras.tau_sq_GHS(k) = ghs_list[k].tau_sq;
          if (idx_k.n_elem == 0) ++empty_precision_updates(k);
          // Restore the original covariance ridge once det(precision) > 1e150.
          arma::mat precision = paras.Sigma_inv.slice(k);
          arma::mat covariance;
          if (ghs_diag_rate == 0.0 && !std::isfinite(ghs_diag_upper) &&
              cortree::regularize_precision(precision, covariance, hyper.err_precision)) {
            paras.Sigma_inv.slice(k) = precision;
            // The next GHS sweep must start from this same inverse pair.
            ghs_list[k].Omega = precision;
            ghs_list[k].Sigma = covariance;
            ++covariance_regularizations(k);
          }
        } else {
          // Original freeze threshold: leave correlated covariance unchanged.
          ++covariance_updates_skipped(k);
        }
      }
      // Legacy flat diagonals retain empty-component precision; either proper
      // diagonal prior refreshes it from the GHS transitions above.
      // Proper independent variance priors can still be refreshed when n_k = 0.
      if (tree.idx_ind.n_elem > 0) {
        arma::mat residual = paras.phi(idx_k, tree.idx_ind);
        residual.each_row() -= paras.mu(tree.idx_ind, arma::uvec{k}).t();
        paras.sigma2_vec.col(k) = cortree::independent_variance(
          residual, tree.ind_layer_idx, hyper.c_sigma2_vec);
      }
    }
    if (!all_ind && ghs_scale_hierarchy)
      precision_scale_rate = cortree::ghs_common_scale_rate(
        precision_scales, ghs_scale_shape, ghs_scale_rate_shape, ghs_scale_rate_rate);
  };

  void update_Z(){ 
    double gamma_det = 1.0;
    if(z_det_gamma_linear && gibbs_control.burnin > 1 && iter < gibbs_control.burnin){
      double frac = static_cast<double>(iter) / static_cast<double>(gibbs_control.burnin - 1);
      gamma_det = z_det_gamma0 + (1.0 - z_det_gamma0) * frac;
      gamma_det = std::min(1.0, std::max(0.0, gamma_det));
    }
    arma::vec log_pi = arma::log(paras.pi);
    if(z_mode == 1){
      // v1-like assignment update (ablation): omit log-det constants and finite guards
      for(arma::uword i=0; i<dat.n; i++){
        arma::vec loglik_cor_i(dat.n_clus); loglik_cor_i.zeros();
        arma::vec loglik_ind_i(dat.n_clus);
        arma::vec phi_sub;
        if(!all_ind){
          phi_sub = arma::trans(paras.phi(arma::uvec{i},tree.idx_cor));
        }
        arma::vec phi_ind = arma::trans(paras.phi(arma::uvec{i},tree.idx_ind));
        for(arma::uword k=0; k<dat.n_clus; k++){
          if(!all_ind){
            arma::vec mu_sub = paras.mu(tree.idx_cor, arma::uvec{k});
            arma::vec diff_sub = phi_sub - mu_sub;
            loglik_cor_i(k) = -0.5*arma::dot(diff_sub, paras.Sigma_inv.slice(k) * diff_sub);
          }
          arma::vec mu_ind = paras.mu(tree.idx_ind, arma::uvec{k});
          arma::vec diff_ind = phi_ind - mu_ind;
          loglik_ind_i(k) = -0.5*arma::dot(diff_ind, diff_ind/paras.sigma2_vec.col(k));
        }
        arma::vec loglik_i = loglik_ind_i + loglik_cor_i;
        arma::vec log_prob = loglik_i + log_pi;
        log_prob -= arma::max(log_prob);
        log_prob = arma::exp(log_prob);
        log_prob /= arma::sum(log_prob);
        arma::vec CDF = arma::cumsum(log_prob);
        double u = arma::randu();
        if(max(CDF) < u){
          Rcout<<"---Error: update_Z::CDF="<<CDF.t()<<std::endl;
        }
        paras.Z(i) = arma::as_scalar(arma::find(CDF >= u, 1));
      }
    }else{
      // current robust assignment update
      arma::vec loglik_cor_const(dat.n_clus, arma::fill::zeros);
      arma::uvec cor_valid(dat.n_clus, arma::fill::ones);
      arma::vec log_det_ind(dat.n_clus, arma::fill::zeros);
      arma::uvec ind_valid(dat.n_clus, arma::fill::ones);
      for(arma::uword k=0; k<dat.n_clus; k++){
        if(!all_ind){
          double log_det_val = 0.0;
          double sign_det = 0.0;
          arma::log_det(log_det_val, sign_det, paras.Sigma_inv.slice(k));
          if(sign_det <= 0.0 || !arma::is_finite(log_det_val)){
            cor_valid(k) = 0;
          }else{
            loglik_cor_const(k) = 0.5 * gamma_det * log_det_val;
          }
        }
        const arma::subview_col<double> sigma2_k = paras.sigma2_vec.col(k);
        if(arma::any(sigma2_k <= 0.0) || !sigma2_k.is_finite()){
          ind_valid(k) = 0;
        }else{
          log_det_ind(k) = arma::sum(arma::log(sigma2_k));
        }
      }

      for(arma::uword i=0; i<dat.n; i++){
        arma::vec loglik_cor_i(dat.n_clus); loglik_cor_i.zeros();
        arma::vec loglik_ind_i(dat.n_clus);
        arma::vec phi_sub;
        if(!all_ind){
          phi_sub = arma::trans(paras.phi(arma::uvec{i},tree.idx_cor));
        }
        arma::vec phi_ind = arma::trans(paras.phi(arma::uvec{i},tree.idx_ind));
        for(arma::uword k=0; k<dat.n_clus; k++){
          if(!all_ind){
            if(cor_valid(k) == 0){
              loglik_cor_i(k) = -arma::datum::inf;
            }else{
              arma::vec mu_sub = paras.mu(tree.idx_cor, arma::uvec{k});
              arma::vec diff_sub = phi_sub - mu_sub;
              loglik_cor_i(k) = -0.5 * arma::dot(diff_sub, paras.Sigma_inv.slice(k) * diff_sub) + loglik_cor_const(k);
            }
          }
          arma::vec mu_ind = paras.mu(tree.idx_ind, arma::uvec{k});
          arma::vec diff_ind = phi_ind - mu_ind;
          const arma::subview_col<double> sigma2_k = paras.sigma2_vec.col(k);
          if(ind_valid(k) == 0){
            loglik_ind_i(k) = -arma::datum::inf;
          }else{
            double quad_ind = arma::dot(diff_ind, diff_ind / sigma2_k);
            loglik_ind_i(k) = -0.5 * (quad_ind + gamma_det * log_det_ind(k));
          }
        }
        arma::vec loglik_i = loglik_ind_i + loglik_cor_i;
        arma::vec log_prob = loglik_i + log_pi;
        arma::uvec finite_idx = arma::find_finite(log_prob);
        if(finite_idx.n_elem == 0){
          log_prob.fill(1.0/static_cast<double>(dat.n_clus));
        }else{
          arma::vec prob = arma::zeros<arma::vec>(dat.n_clus);
          double max_finite = log_prob(finite_idx).max();
          prob(finite_idx) = arma::exp(log_prob(finite_idx) - max_finite);
          double prob_sum = arma::sum(prob);
          if(!arma::is_finite(prob_sum) || prob_sum <= 0.0){
            log_prob.fill(1.0/static_cast<double>(dat.n_clus));
          }else{
            log_prob = prob / prob_sum;
          }
        }
        arma::vec CDF = arma::cumsum(log_prob);
        double u = arma::randu();
        if(max(CDF) < u){
          Rcout<<"---Error: update_Z::CDF="<<CDF.t()<<std::endl;
        }
        paras.Z(i) = arma::as_scalar(arma::find(CDF >= u, 1));
      }
    }

    // release the covariance status
    for(int k=0; k<dat.n_clus; k++){
      arma::uvec idx_k = arma::find(paras.Z == k);
      if(idx_k.n_elem >= 30){
        paras.Z_cor_status(k) = 1;
      }
    }
    
    
  };
  
  void update_pi(){
    arma::uvec counts(dat.n_clus, arma::fill::zeros);
    for (arma::uword i = 0; i < paras.Z.n_elem; ++i) ++counts(paras.Z(i));
    paras.pi = cortree::stick_weights(counts, paras.alpha);
  };

  void update_loglike(){
    const arma::mat parent_count = tree.count.cols(0, tree.total_parents - 1);
    paras.loglike = 0.0;
    for (arma::uword i = 0; i < parent_count.n_rows; ++i) {
      for (arma::uword j = 0; j < parent_count.n_cols; ++j) {
        double n = parent_count(i,j);
        paras.loglike += cortree::binomial_log_kernel(n, tree.kappa(i,j) + n / 2.0, paras.phi(i,j));
      }
    }

    if(save_cluster_cor_trace){
      // update cluster correlation
      arma::umat indicator_matrix(dat.n, dat.n, arma::fill::zeros);
      for (arma::uword i = 0; i < dat.n; ++i) {
        indicator_matrix.col(i) = (paras.Z == paras.Z(i));
      }
      paras.cluster_cor = indicator_matrix;
    }
  }
  
  void run_gibbs(){
    for(iter=0; iter<gibbs_control.total_iter; iter++){
      Rcpp::checkUserInterrupt();
      update_omega();
      update_phi();
      update_mu();

      if(all_ind){
        update_Sigma();
      }else if(iter % cov_interval ==0){
        update_Sigma();
      }
      
      
      update_Z();
      update_pi();
      
      
      update_loglike();
      save_gibbs_sample();

      // Old wrong code: `total_iter < 10` makes `ten_percent = 0`,
      // and `iter % ten_percent` triggers divide-by-zero/undefined behavior.
      // int ten_percent = gibbs_control.total_iter/10;
      // if(iter % ten_percent == 0){
      int ten_percent = std::max(1, gibbs_control.total_iter / 10);
      if(iter % ten_percent == 0){
        Rcpp::Rcout << "iter: " << iter << " loglike: " << paras.loglike << std::endl;
        Rcout<<"---update_pi::pi="<<paras.pi.t()<<std::endl;
        // Rcout<<"status of covariance ="<<paras.Z_cor_status.t()<<std::endl;
      }
    }
  };
  
  void save_gibbs_sample(){
    paras_sample.pi_sample.col(iter) = paras.pi;
    // Old wrong code: this skips `iter == burnin`, so only
    // (total_iter - burnin - 1) draws are stored even though memory is
    // allocated for (total_iter - burnin) draws.
    // if(iter > gibbs_control.burnin){
    //   int idx = iter - gibbs_control.burnin - 1;
    if(iter >= gibbs_control.burnin){
      int idx = iter - gibbs_control.burnin;
      paras_sample.mu_sample.slice(idx) = paras.mu;
      if(!all_ind){
        paras_sample.Sigma_inv_sample[idx] = paras.Sigma_inv;
      }
      if(save_phi_trace){
        paras_sample.phi_sample.slice(idx) = paras.phi;
      }
      paras_sample.Z_sample.col(idx) = paras.Z;
      if(save_cluster_cor_trace){
        paras_sample.cluster_cor.slice(idx) = paras.cluster_cor;
      }
      paras_sample.tausq_GHS.col(idx) = paras.tau_sq_GHS;
      if (!all_ind && ghs_scale_hierarchy) {
        paras_sample.scale_t.col(idx) = precision_scales;
        paras_sample.scale_b(idx) = precision_scale_rate;
      }
      paras_sample.sigma2_vec_sample.slice(idx) = paras.sigma2_vec;
    }
    paras_sample.loglik(iter) = paras.loglike;
  };

  Rcpp::List get_gibbs_sample(){
    SEXP phi_output = save_phi_trace ? Rcpp::wrap(paras_sample.phi_sample) : R_NilValue;
    SEXP cluster_cor_output = save_cluster_cor_trace ? Rcpp::wrap(paras_sample.cluster_cor) : R_NilValue;
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("mu") = paras_sample.mu_sample,
                              Rcpp::Named("phi") = phi_output,
                              Rcpp::Named("cluster_cor") = cluster_cor_output,
                              Rcpp::Named("Sigma_inv") = paras_sample.Sigma_inv_sample,
                              Rcpp::Named("sigma2_vec") = paras_sample.sigma2_vec_sample,
                              Rcpp::Named("tausq_GHS") = paras_sample.tausq_GHS,
                              Rcpp::Named("pi") = paras_sample.pi_sample,
                              Rcpp::Named("Z") = paras_sample.Z_sample,
                              Rcpp::Named("loglik") = paras_sample.loglik,
                              Rcpp::Named("ghs_prior") = Rcpp::List::create(
                                Rcpp::Named("family") = ghs_scale_hierarchy ? "hierarchical_determinant_trace" :
                                  ghs_det_df > 0.0 ? (std::isfinite(ghs_diag_upper) ? "bounded_determinant" : "determinant_trace") :
                                  ghs_jmlr_lambda > 0.0 ? "jmlr_fixed_scale" : "random_global_scale",
                                Rcpp::Named("det_df") = ghs_det_df,
                                Rcpp::Named("determinant_power") = ghs_det_df / 2.0,
                                Rcpp::Named("jmlr_lambda") = ghs_jmlr_lambda,
                                Rcpp::Named("global_tau_fixed") = ghs_jmlr_lambda > 0.0,
                                Rcpp::Named("global_scale_prior") = ghs_jmlr_lambda > 0.0 ? "fixed" : "unit_half_cauchy",
                                Rcpp::Named("global_scale") = ghs_jmlr_lambda > 0.0 ? 1.0 / ghs_jmlr_lambda : NA_REAL,
                                Rcpp::Named("local_scale_prior") = "unit_half_cauchy",
                                Rcpp::Named("diagonal") = ghs_det_df > 0.0 ?
                                  (std::isfinite(ghs_diag_upper) ? "bounded_determinant_tilt" : "exponential_determinant_tilt") :
                                  std::isfinite(ghs_diag_upper) ? "uniform" :
                                  (ghs_diag_rate > 0.0 ? "exponential" : "flat_legacy"),
                                Rcpp::Named("diag_rate") = ghs_diag_rate,
                                Rcpp::Named("diag_upper") = ghs_diag_upper,
                                Rcpp::Named("proper") = ghs_diag_rate > 0.0 || std::isfinite(ghs_diag_upper),
                                Rcpp::Named("active") = !all_ind,
                                Rcpp::Named("normalization") = "joint_over_SPD_and_horseshoe_scales",
                                Rcpp::Named("warm_start") = warm_start),
                              Rcpp::Named("empty_precision_updates") = empty_precision_updates,
                              Rcpp::Named("ghs_uniform") = Rcpp::List::create(
                                Rcpp::Named("block_updates") = uniform_block_updates,
                                Rcpp::Named("ess_evaluations") = uniform_ess_evaluations,
                                Rcpp::Named("max_ess_evaluations") = uniform_max_ess_evaluations),
                              Rcpp::Named("covariance_safeguards") = Rcpp::List::create(
                                Rcpp::Named("enabled") = !all_ind && ghs_diag_rate == 0.0 && !std::isfinite(ghs_diag_upper),
                                Rcpp::Named("regularizations") = covariance_regularizations,
                                Rcpp::Named("updates_skipped") = covariance_updates_skipped,
                                Rcpp::Named("epsilon") = hyper.err_precision,
                                Rcpp::Named("regularize_logdet_threshold") = std::log(1e150),
                                Rcpp::Named("freeze_logdet_threshold") = std::log(1e200)));
    if (ghs_scale_hierarchy) {
      Rcpp::List prior = output["ghs_prior"];
      prior["scale_hierarchy"] = true;
      prior["precision_parameterization"] = "Omega=t*Q";
      prior["diag_rate_applies_to"] = "Q";
      prior["horseshoe_scales_apply_to"] = "Q";
      prior["template_diag_rate"] = ghs_diag_rate;
      prior["template_det_df"] = ghs_det_df;
      output["ghs_prior"] = prior;
      output["ghs_scale"] = Rcpp::List::create(
        Rcpp::Named("active") = !all_ind,
        Rcpp::Named("t") = !all_ind ? Rcpp::wrap(paras_sample.scale_t) : R_NilValue,
        Rcpp::Named("b") = !all_ind ? Rcpp::wrap(paras_sample.scale_b) : R_NilValue,
        Rcpp::Named("shape") = ghs_scale_shape,
        Rcpp::Named("rate_shape") = ghs_scale_rate_shape,
        Rcpp::Named("rate_rate") = ghs_scale_rate_rate,
        Rcpp::Named("initial_t") = 1.0,
        Rcpp::Named("initial_b") = ghs_scale_shape > 1.0 ? ghs_scale_shape - 1.0 : ghs_scale_shape,
        Rcpp::Named("trace_scope") = "post_burnin",
        Rcpp::Named("first_saved_iteration") = gibbs_control.burnin + 1,
        Rcpp::Named("update_schedule") = "cov_interval",
        Rcpp::Named("update_interval") = cov_interval,
        Rcpp::Named("Gamma_parameterization") = "shape_rate",
        Rcpp::Named("precision_trace") = "Sigma_inv=t*Q",
        Rcpp::Named("global_scale_trace") = "tausq_GHS=template_Q_global_variance");
    }
    return output;
  };
  
  
};

// The explicit R signature preserves Inf and required init_Z when regenerating exports.
// [[Rcpp::export(signature = {X, n_clus, tree_depth, cutoff_layer, total_iter, burnin, warm_start = 0L, init_Z, c_sigma2_vec = 1.0, sigma_mu2 = 1.0, all_ind = FALSE, cov_interval = 1L, save_phi_trace = FALSE, save_cluster_cor_trace = FALSE, z_det_gamma0 = 1.0, z_mode = 0L, z_det_gamma_linear = FALSE, ghs_diag_rate = 0.0, ghs_diag_upper = Inf, ghs_jmlr_lambda = 0.0, ghs_det_df = 0.0, ghs_scale_hierarchy = FALSE, ghs_scale_shape = 3.0, ghs_scale_rate_shape = 2.0, ghs_scale_rate_rate = 1.0})]]
Rcpp::List CorTree_sampler(arma::mat X, 
                      int n_clus, int tree_depth, int cutoff_layer, 
                      int total_iter, int burnin, int warm_start=0,
                      arma::uvec init_Z = arma::zeros<arma::uvec>(1),
                      double c_sigma2_vec = 1.0, 
                      double sigma_mu2=1.0,
                      bool all_ind = false,
                      int cov_interval = 1,
                      bool save_phi_trace = false,
                      bool save_cluster_cor_trace = false,
                      double z_det_gamma0 = 1.0,
                      int z_mode = 0,
                      bool z_det_gamma_linear = false,
                      double ghs_diag_rate = 0.0,
                      double ghs_diag_upper = R_PosInf,
                      double ghs_jmlr_lambda = 0.0,
                      double ghs_det_df = 0.0,
                      bool ghs_scale_hierarchy = false,
                      double ghs_scale_shape = 3.0,
                      double ghs_scale_rate_shape = 2.0,
                      double ghs_scale_rate_rate = 1.0){
  arma::wall_clock timer;
  timer.tic();
  cortree::validate_sampler(X, n_clus, total_iter, burnin, warm_start,
                            cov_interval, c_sigma2_vec, sigma_mu2, init_Z);
  cortree::validate_ghs_diag_rate(ghs_diag_rate);
  cortree::validate_ghs_jmlr_lambda(ghs_jmlr_lambda, ghs_diag_rate, ghs_diag_upper);
  cortree::validate_ghs_det_df(ghs_det_df, ghs_diag_upper, ghs_diag_rate, ghs_jmlr_lambda);
  cortree::validate_ghs_scale_hierarchy(ghs_scale_hierarchy, ghs_scale_shape,
    ghs_scale_rate_shape, ghs_scale_rate_rate, ghs_diag_rate, ghs_diag_upper,
    ghs_jmlr_lambda, warm_start);
  cortree::validate_ghs_diag_upper(ghs_diag_upper, ghs_diag_rate, warm_start, all_ind);
  cortree::validate_depth(tree_depth, cutoff_layer, all_ind);
  if (z_mode != 0 && z_mode != 1) Rcpp::stop("z_mode must be 0 or 1.");
  if (!std::isfinite(z_det_gamma0) || z_det_gamma0 < 0 || z_det_gamma0 > 1)
    Rcpp::stop("z_det_gamma0 must be in [0,1].");
  if (z_mode == 1) Rcpp::warning("z_mode=1 is a legacy ablation, not a posterior sampler.");
  CorTree model;

  if(init_Z.n_elem == 1 && X.n_rows > 1){
    init_Z = arma::randi<arma::uvec>(X.n_rows, arma::distr_param(0, n_clus-1));
  }
  
  model.all_ind = all_ind;
  model.save_phi_trace = save_phi_trace;
  model.save_cluster_cor_trace = save_cluster_cor_trace;
  model.warm_start = warm_start;
  model.cov_interval = cov_interval;
  model.ghs_jmlr_lambda = ghs_jmlr_lambda;
  model.ghs_det_df = ghs_det_df;
  model.ghs_scale_hierarchy = ghs_scale_hierarchy;
  model.ghs_scale_shape = ghs_scale_shape;
  model.ghs_scale_rate_shape = ghs_scale_rate_shape;
  model.ghs_scale_rate_rate = ghs_scale_rate_rate;
  model.ghs_diag_rate = ghs_jmlr_lambda > 0.0 ? ghs_jmlr_lambda / 2.0 : ghs_diag_rate;
  model.ghs_diag_upper = ghs_diag_upper;
  model.z_det_gamma0 = z_det_gamma0;
  model.z_mode = z_mode;
  model.z_det_gamma_linear = z_det_gamma_linear;
  model.load_data(X, n_clus);
  Rcout << "Data loaded" << std::endl;
  model.set_hyperparameter(c_sigma2_vec, sigma_mu2);
  Rcout << "Hyperparameter set" << std::endl;
  model.set_gibbs_control(total_iter, burnin);
  Rcout << "Gibbs control set" << std::endl;
  model.vectorize_tree(tree_depth, cutoff_layer);
  Rcout << "Tree vectorized" << std::endl;
  model.initialize_parameters(init_Z);
  Rcout << "Parameters initialized" << std::endl;
  model.initialize_mcmc_sample();
  Rcout << "MCMC sample initialized" << std::endl;
  model.run_gibbs();
  Rcout << "Gibbs sampler finished" << std::endl;
  double elapsed = timer.toc();
  
  Rcpp::List output;
  model.set_test_output();
  output = Rcpp::List::create(Rcpp::Named("mcmc") = model.get_gibbs_sample(),
                              Rcpp::Named("test_output") = model.test_output,
                              Rcpp::Named("elapsed") = elapsed);
  
  return output;
}
