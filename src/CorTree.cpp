#include "PolyaGamma.h" // to sample polya-gamma random variable
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>



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
  } paras_sample;
  
public:
  
  int iter=0;
  std::vector<GHS_params> ghs_list;
  Rcpp::List test_output;
  bool all_ind = false;
  int warm_start = 0;
  int cov_interval = 1;
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
  
  void initialize_mcmc_sample(){
    paras_sample.mu_sample = arma::zeros<arma::cube>(tree.total_parents, dat.n_clus, gibbs_control.mcmc_sample);
    paras_sample.phi_sample = arma::zeros<arma::cube>(dat.n, tree.total_parents, gibbs_control.mcmc_sample);
    if(!all_ind){
      paras_sample.Sigma_inv_sample = Rcpp::List(gibbs_control.mcmc_sample);
    }
    paras_sample.pi_sample = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.mcmc_sample);
    paras_sample.loglik = arma::zeros<arma::vec>(gibbs_control.total_iter);
    paras_sample.Z_sample = arma::zeros<arma::umat>(dat.n, gibbs_control.mcmc_sample);
    
    paras_sample.cluster_cor = arma::zeros<arma::ucube>(dat.n, dat.n, gibbs_control.mcmc_sample);
    paras_sample.tausq_GHS = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.mcmc_sample);
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
      tree.idx_ind = arma::regspace<arma::uvec>(tree.L, tree.total_parents-1);
    }
    
    
    // ser precision for Sigma_inv
    hyper.err_precision = pow(10.0, -15.0);
    Rcout<<"err_precision="<<hyper.err_precision<<std::endl;
    
    // initialize the tree
    tree.count = arma::zeros<arma::mat>(dat.X.n_rows, tree.total_nodes);
    tree.kappa = arma::zeros<arma::mat>(dat.X.n_rows, tree.total_parents);
    tree.interval_left = arma::zeros<arma::uvec>(tree.total_nodes);
    tree.interval_right = arma::zeros<arma::uvec>(tree.total_nodes);
    
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
      if(tree.interval_left(i) <= tree.interval_right(left_node)-1){
        left_idx = arma::regspace<arma::uvec>(tree.interval_left(i), tree.interval_right(left_node)-1);
        tree.count.col(left_node) = arma::sum(dat.X.cols(left_idx), 1);
      }else{
        tree.count.col(left_node) = arma::zeros<arma::vec>(dat.X.n_rows);
      }
      
      
      if(tree.interval_left(right_node) <= tree.interval_right(i)-1){
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
    if(!all_ind){
      tree.ind_layer_idx = floor( log2(arma::regspace<arma::vec>(tree.L, tree.total_parents - 1)) );
    }else{
      tree.ind_layer_idx = floor( log2(2+arma::regspace<arma::vec>(0, tree.total_parents - 1)) );
    }
    
  }; 
  
  void set_test_output(){
    test_output = Rcpp::List::create(Rcpp::Named("count") = tree.count,
                                     Rcpp::Named("sigma_vec") = paras.sigma2_vec,
                                     Rcpp::Named("kappa") = tree.kappa,
                                     Rcpp::Named("interval_left") = tree.interval_left,
                                     Rcpp::Named("interval_right") = tree.interval_right);
  }
  
  
  void initialize_parameters(arma::uvec init_Z){
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
    }
    
    paras.tau_sq_GHS = arma::zeros<arma::vec>(dat.n_clus);
    
  };
  
  
  
  
  void update_omega(){
    
    // #pragma omp parallel for
    for(int i = 0; i < dat.X.n_rows; i++){
      arma::rowvec b_i = tree.count.row(i).cols(0, tree.total_parents-1);
      arma::rowvec c_i = paras.phi.row(paras.Z(i));
      // arma::rowvec omega_i = arma_pgdraw(b_i, c_i);
      NumericVector omega_i_vec = rcpp_pgdraw(Rcpp::NumericVector(b_i.begin(), b_i.end()), 
                                              Rcpp::NumericVector(c_i.begin(), c_i.end()));
      paras.omega.row(i) = arma::rowvec(omega_i_vec.begin(), omega_i_vec.size());
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
        arma::mat temp;
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
      
      // to avoid numerical issue, phi_i cannot be larger than 10
      phi_i.clamp( -7.0, 7.0);
      
      paras.phi.row(i) = phi_i.t();
    }
    
  };
  
  void update_mu(){
    // #pragma omp parallel for
    for(arma::uword k=0; k<dat.n_clus; k++){
      arma::uvec idx_k = arma::find(paras.Z == k);
      if(idx_k.n_elem > 0){
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
    
    // arma::mat cholDecomp;
    // bool pd_check = arma::chol(cholDecomp, params.Omega);
    arma::mat Omega_old = params.Omega;
    // bool S_pd_check = arma::chol(cholDecomp, S);
    
    // Indices for columns
    arma::mat ind_all(p - 1, p, arma::fill::zeros);
    for (int i = 0; i < p; ++i) {
      if (i == 0) {
        ind_all.col(i) = arma::regspace(1, p - 1);
      } else if (i == p - 1) {
        ind_all.col(i) = arma::regspace(0, p - 2);
      } else {
        arma::uvec inds = arma::regspace<arma::uvec>(0, p - 1);
        inds.shed_row(i);
        ind_all.col(i) = arma::conv_to<arma::vec>::from(inds);
      }
    }
    
    
    // Sample Sigma and Omega = inv(Sigma)
    for (arma::uword i = 0; i < p; ++i) {
      arma::uvec ind = arma::conv_to<arma::uvec>::from(ind_all.col(i));
      
      // Extract submatrices using proper matrix subsetting
      arma::mat Sigma_11 = params.Sigma(ind, ind);             // Sigma_11
      arma::vec sigma_12 = params.Sigma(ind, arma::uvec{i});    // Sigma_12
      double sigma_22 = params.Sigma(i,i);  // Sigma_22 as submatrix
      
      arma::vec s_21 = S(ind, arma::uvec{i});            // s_21
      double s_22 = S(i,i);  // s_22 as submatrix
      
      arma::vec lambda_sq_12 = params.Lambda_sq(ind, arma::uvec{i});    // Lambda_sq_12
      arma::vec nu_12 = params.Nu(ind, arma::uvec{i});                  // Nu_12
      
      // Sample gamma and beta using arma::randg
      double gamma = arma::randg<double>(arma::distr_param(n / 2.0 + 1, 2.0 / s_22));  // Convert s_22 to scalar
      arma::mat inv_Omega_11 = Sigma_11 - sigma_12 * sigma_12.t() / sigma_22;  // Convert sigma_22 to scalar
      
      arma::mat inv_C = s_22 * inv_Omega_11 + arma::diagmat(1.0 / (lambda_sq_12 * params.tau_sq));
      arma::mat inv_C_chol = arma::chol(inv_C);
      arma::vec mu_i = -arma::solve(inv_C, s_21);
      arma::vec beta = mu_i + arma::solve(inv_C_chol, arma::randn(p - 1));
      arma::vec omega_12 = beta;
      double omega_22 = gamma + arma::as_scalar(beta.t() * inv_Omega_11 * beta);
      
      // Update Lambda_sq and Nu using arma::randg
      arma::vec rate = arma::square(omega_12) / (2.0 * params.tau_sq) + 1.0 / nu_12;
      arma::mat inv_lambda_sq_12 = arma::randg<arma::vec>(rate.n_elem, arma::distr_param(1.0,1.0)); 
      inv_lambda_sq_12 /= rate;
      lambda_sq_12 = 1.0 / inv_lambda_sq_12;
      nu_12 = arma::randg<arma::vec>(nu_12.n_elem, arma::distr_param(1.0,1.0));
      nu_12 /= lambda_sq_12;
      nu_12 = 1.0 / nu_12;
      
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
      
      // Update tau_sq and xi
      arma::vec omega_vector = params.Omega(arma::trimatl_ind(size(params.Omega), -1));
      arma::vec lambda_sq_vector = params.Lambda_sq(arma::trimatl_ind(size(params.Lambda_sq), -1));
      double rate_tau_sq = 1.0 / params.xi + arma::sum(arma::square(omega_vector) / (2.0 * lambda_sq_vector));
      params.tau_sq = 1.0 / arma::randg<double>(arma::distr_param((p * (p - 1) / 2 + 1) / 2.0, 1.0 / rate_tau_sq));
      params.xi = 1.0 / arma::randg<double>(arma::distr_param(1.0, 1.0 / (1.0 + 1.0 / params.tau_sq)));
    }
    
    
  }
  
  
  void update_Sigma(){
    // consider a warm start using independent state
    for(arma::uword k=0; k<dat.n_clus; k++){
      arma::uvec idx_k = arma::find(paras.Z == k);
      if(idx_k.n_elem > 0){
        if(!all_ind){
          if(iter < warm_start || paras.Z_cor_status(k) == 0){
            Rcout<<"---update_Sigma::warm start for "<<k<<"-th cluster, Z_cor_status="<<paras.Z_cor_status(k)<<std::endl;
            
            // // update sigma_vec
            arma::vec layer_idx = floor( log2(2+arma::regspace<arma::vec>(0, tree.L - 1)) );
            // arma::vec a_vec = layer_idx* hyper.c_sigma2_vec + idx_k.n_elem/2;
            // arma::mat phi_res = paras.phi(idx_k,tree.idx_cor);
            // phi_res.each_row() -= paras.mu(tree.idx_cor, arma::uvec{k}).t();
            // arma::vec b_vec = layer_idx + arma::trans(arma::sum(phi_res%phi_res, 0)/2);  
            
            // // wrong
            // arma::vec gamma1(a_vec.n_elem);
            // arma::vec gamma2(b_vec.n_elem);
            // for (arma::uword i = 0; i < a_vec.n_elem; ++i) {
            //   gamma1(i) = arma::randg(arma::distr_param(a_vec(i), 1.0));  // generate gamma random values
            //   gamma2(i) = arma::randg(arma::distr_param(b_vec(i), 1.0));
            // }
            // arma::vec sigma2_k_vec = (gamma1 + gamma2)/gamma1;
            // paras.Sigma_inv.slice(k) = arma::diagmat(sigma2_k_vec);
            
            // correct
            double a_double = hyper.c_sigma2_vec + idx_k.n_elem/2;
            arma::mat phi_res = paras.phi(idx_k,tree.idx_cor);
            phi_res.each_row() -= paras.mu(tree.idx_cor, arma::uvec{k}).t();
            arma::vec phi_mse = arma::trans(arma::sum(phi_res%phi_res, 0)/2);
            // arma::vec b_vec = tree.ind_layer_idx + phi_mse;  
            arma::vec b_vec = 1.0/layer_idx + phi_mse;  
            arma::vec gamma(b_vec.n_elem);
            for (arma::uword i = 0; i < b_vec.n_elem; ++i) {
              // gamma(i) = arma::randg(arma::distr_param(a_vec(i), 1.0/b_vec(i)));  // generate gamma random values
              gamma(i) = arma::randg(arma::distr_param(a_double, 1.0/b_vec(i)));  // generate gamma random values
            }
            
            paras.Sigma_inv.slice(k) = arma::diagmat(gamma);
            
            
          }else if(arma::det(paras.Sigma_inv.slice(k)) < 1e200){
            
            arma::mat phi_centered = paras.phi(idx_k, tree.idx_cor);
            phi_centered.each_row() -= paras.mu(tree.idx_cor, arma::uvec{k}).t();
            // arma::rowvec phi_mean = arma::mean(phi_centered, 0);
            // phi_centered.each_row() -= phi_mean;
            arma::mat phi_centered_Cov = phi_centered.t() * phi_centered;
            GHS_oneSample(phi_centered_Cov, dat.n, ghs_list[k]);
            paras.tau_sq_GHS(k) = ghs_list[k].tau_sq;
            // check if Omega is symmetric positive definite
            
            paras.Sigma_inv.slice(k) = ghs_list[k].Omega;
            
            // avoids singular updates
            if(arma::det(paras.Sigma_inv.slice(k)) > 1e150){
              // Eigen-value regularization
              Rcout<<"Eigen-value regularization: adding a small pertubation "<<hyper.err_precision<<" to "<<k<<"-th covariance"<<std::endl;
              // Eigen decomposition of Sigma
              arma::vec eigvals;
              arma::mat eigvecs;
              eig_sym(eigvals, eigvecs, paras.Sigma_inv.slice(k));
              eigvals = 1.0/eigvals + hyper.err_precision;  // Add epsilon to each eigenvalue
              arma::vec eigvals_inv = 1.0 / eigvals;
              paras.Sigma_inv.slice(k) = eigvecs * diagmat(eigvals_inv) * eigvecs.t();
              Rcout<<" arma::det(paras.Sigma_inv.slice(k)) ="<<arma::det(paras.Sigma_inv.slice(k)) <<std::endl;
            }
            // paras.Sigma_inv.slice(k).diag() += hyper.err_precision*arma::ones(tree.L); // to prevent singularity
            
            // double prop_k = idx_k.n_elem/dat.n;
            // if( idx_k.n_elem < 30){
            //   Rcout<<"make cov ind: counts of idx_"<<k<<" is "<<idx_k.n_elem<<"; idx_k.n_elem="<<idx_k.n_elem <<std::endl;
            //   paras.Z_cor_status(k) = 0;
            // }
            
            if(iter % cov_interval == 0){
              Rcout<<"---update_Sigma::det of Sigma_inv["<<k<<"]="<<arma::det(paras.Sigma_inv.slice(k))<<std::endl;
            }
          }
          
          
        }
        
        // update sigma_vec
        // arma::vec a_vec = tree.ind_layer_idx * hyper.c_sigma2_vec + idx_k.n_elem/2;
        double a_double = hyper.c_sigma2_vec + idx_k.n_elem/2;
        arma::mat phi_res = paras.phi(idx_k,tree.idx_ind);
        phi_res.each_row() -= paras.mu(tree.idx_ind, arma::uvec{k}).t();
        arma::vec phi_mse = arma::trans(arma::sum(phi_res%phi_res, 0)/2);
        // arma::vec b_vec = tree.ind_layer_idx + phi_mse;  
        arma::vec b_vec = 1.0/tree.ind_layer_idx + phi_mse;  
        
        // // wrong
        // arma::vec gamma1(a_vec.n_elem);
        // arma::vec gamma2(b_vec.n_elem);
        // for (arma::uword i = 0; i < a_vec.n_elem; ++i) {
        //   gamma1(i) = arma::randg(arma::distr_param(a_vec(i), 1.0));  // generate gamma random values
        //   gamma2(i) = arma::randg(arma::distr_param(b_vec(i), 1.0));
        // }
        // paras.sigma2_vec.col(k) = (gamma1 + gamma2)/gamma1;
        
        // correct: inverse gamma
        arma::vec gamma(b_vec.n_elem);
        for (arma::uword i = 0; i < b_vec.n_elem; ++i) {
          // gamma(i) = arma::randg(arma::distr_param(a_vec(i), 1.0/b_vec(i)));  // generate gamma random values
          gamma(i) = arma::randg(arma::distr_param(a_double, 1.0/b_vec(i)));  // generate gamma random values
        }
        paras.sigma2_vec.col(k) = 1.0/gamma;
        
        // // correct: half-cauchy
        // arma::vec gamma1(b_vec.n_elem);
        // arma::vec gamma_a(b_vec.n_elem);
        // arma::vec inv_a_k = paras.inv_a_vec.col(k);
        // for (arma::uword i = 0; i < b_vec.n_elem; ++i) {
        //   gamma1(i) = arma::randg(arma::distr_param(idx_k.n_elem/2+0.5, 
        //                           1.0/(phi_mse(i) + inv_a_k(i)) ));  // generate gamma random values
        //   gamma_a(i) = arma::randg(arma::distr_param(1.0,1.0/(1.0/gamma1(i) + tree.ind_layer_idx(i))));  // generate gamma random values
        // }
        // paras.sigma2_vec.col(k) = 1.0/gamma1;
        // paras.inv_a_vec = 1.0/gamma_a;
        
        
        
      }
      
    }
    
  };// Graphical Horseshoe prior
  
  void update_Z(){ 
    // #pragma omp parallel for
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
      
      arma::vec log_prob = loglik_i + log(paras.pi);
      
      log_prob -= arma::max(log_prob);
      log_prob = arma::exp(log_prob);
      log_prob /= arma::sum(log_prob);
      arma::vec CDF = arma::cumsum(log_prob);
      // Draw a sample from the multinomial distribution using the CDF
      double u = arma::randu();
      if(max(CDF) < u){
        Rcout<<"---Error: update_Z::CDF="<<CDF.t()<<std::endl;
      }
      paras.Z(i) = arma::as_scalar(arma::find(CDF >= u, 1));
      
      
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
    // stick-breaking weights
    arma::uvec Z_count =  arma::hist(paras.Z, arma::regspace<arma::uvec>(0, dat.n_clus-1));
    arma::vec beta_a = 1.0 + arma::conv_to<arma::vec>::from(Z_count);
    arma::vec beta_b = paras.alpha + arma::cumsum(arma::reverse(beta_a));
    beta_b = arma::reverse(beta_b);
    beta_b = beta_b.subvec(1, beta_b.n_elem - 1);
    beta_a = beta_a.subvec(0, beta_a.n_elem - 2);
    
    arma::vec gamma1(beta_a.n_elem);
    arma::vec gamma2(beta_b.n_elem);
    for (arma::uword i = 0; i < beta_a.n_elem; ++i) {
      gamma1(i) = arma::randg(arma::distr_param(beta_a(i), 1.0));  // generate gamma random values
      gamma2(i) = arma::randg(arma::distr_param(beta_b(i), 1.0));
    }
    arma::vec V = gamma1/(gamma2+gamma1);
    paras.pi(0) = V(0);
    for(int k=1; k<dat.n_clus-1; k++){
      paras.pi(k) = V(k) * arma::prod(1 - V.subvec(0, k-1));
    }
    paras.pi(dat.n_clus-1) = 1 - arma::sum(paras.pi.subvec(0, dat.n_clus-2));
  }; // Dirichlet process mixture prior
  
  void update_loglike(){
    arma::mat n_left_child = tree.kappa + tree.count.cols(0,tree.total_parents-1)/2;
    arma::mat V = 1.0/(1.0 + arma::exp(-paras.phi));
    arma::mat loglike_mat = n_left_child % arma::log(V) + (tree.count.cols(0,tree.total_parents-1) - n_left_child) % arma::log(1.0-V);
    paras.loglike = arma::accu(loglike_mat);
    
    // update cluster correlation
    arma::umat indicator_matrix(dat.n, dat.n, arma::fill::zeros);
    for (arma::uword i = 0; i < dat.n; ++i) {
      indicator_matrix.col(i) = (paras.Z == paras.Z(i));
    }
    paras.cluster_cor = indicator_matrix;
  }
  
  void run_gibbs(){
    for(iter=0; iter<gibbs_control.total_iter; iter++){
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
      
      int ten_percent = gibbs_control.total_iter/10;
      if(iter % ten_percent == 0){
        Rcpp::Rcout << "iter: " << iter << " loglike: " << paras.loglike << std::endl;
        Rcout<<"---update_pi::pi="<<paras.pi.t()<<std::endl;
        // Rcout<<"status of covariance ="<<paras.Z_cor_status.t()<<std::endl;
      }
    }
  };
  
  void save_gibbs_sample(){
    if(iter > gibbs_control.burnin){
      int idx = iter - gibbs_control.burnin - 1;
      paras_sample.mu_sample.slice(idx) = paras.mu;
      if(!all_ind){
        paras_sample.Sigma_inv_sample[idx] = paras.Sigma_inv;
      }
      paras_sample.pi_sample.col(idx) = paras.pi;
      paras_sample.phi_sample.slice(idx) = paras.phi;
      paras_sample.Z_sample.col(idx) = paras.Z;
      paras_sample.cluster_cor.slice(idx) = paras.cluster_cor;
      paras_sample.tausq_GHS.col(idx) = paras.tau_sq_GHS;
      paras_sample.sigma2_vec_sample.slice(idx) = paras.sigma2_vec;
    }
    paras_sample.loglik(iter) = paras.loglike;
  };
  
  Rcpp::List get_gibbs_sample(){
    return Rcpp::List::create(Rcpp::Named("mu") = paras_sample.mu_sample,
                              Rcpp::Named("phi") = paras_sample.phi_sample,
                              Rcpp::Named("Sigma_inv") = paras_sample.Sigma_inv_sample,
                              Rcpp::Named("sigma2_vec") = paras_sample.sigma2_vec_sample,
                              Rcpp::Named("cluster_cor") = paras_sample.cluster_cor,
                              Rcpp::Named("tausq_GHS") = paras_sample.tausq_GHS,
                              Rcpp::Named("pi") = paras_sample.pi_sample,
                              Rcpp::Named("Z") = paras_sample.Z_sample,
                              Rcpp::Named("loglik") = paras_sample.loglik);
  };
  
  
};

// [[Rcpp::export]]
Rcpp::List CorTree_sampler(arma::mat X, 
                            int n_clus, int tree_depth, int cutoff_layer, 
                            int total_iter, int burnin, int warm_start=0,
                            Rcpp::Nullable<arma::uvec> init_Z_ = R_NilValue,
                            double c_sigma2_vec = 1.0, 
                            double sigma_mu2=1.0,
                            bool all_ind = false,
                            int cov_interval = 1){
  arma::uvec init_Z = init_Z_.isNull()
  ? arma::zeros<arma::uvec>(1)
    : Rcpp::as<arma::uvec>(init_Z_);
  
  arma::wall_clock timer;
  timer.tic();
  CorTree model;
  
  if(init_Z.n_elem == 1){
    init_Z = arma::randi<arma::uvec>(X.n_rows, arma::distr_param(0, n_clus-1));
  }
  
  model.all_ind = all_ind;
  model.warm_start = warm_start;
  model.cov_interval = cov_interval;
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

// [[Rcpp::export]]
Rcpp::List construct_tree(arma::mat& X, arma::uword tree_depth){
  // set tree structure
  int m = tree_depth;
  int total_nodes = pow(2, m+1) - 1;
  int total_parents = pow(2, m) - 1;
  
  
  
  // initialize the tree
  arma::mat count = arma::zeros<arma::mat>(X.n_rows, total_nodes);
  arma::mat empirical_phi = arma::zeros<arma::mat>(X.n_rows, total_parents);
  
  arma::uvec interval_left = arma::zeros<arma::uvec>(total_nodes);
  arma::uvec interval_right = arma::zeros<arma::uvec>(total_nodes);
  
  // fill in the tree by the count matrix
  int k = X.n_cols;
  count.col(0) = arma::sum(X, 1);
  interval_left(0) = 0; // [left, right) left closed, right open
  interval_right(0) = k;
  
  for(int i=0; i<total_parents; i++){
    
    arma::uword left_node = 2*i+1;
    interval_left(left_node) = interval_left(i);
    interval_right(left_node) = (interval_left(i) + interval_right(i))/2;
    
    arma::uword right_node = 2*i+2;
    interval_left(right_node) = interval_right(left_node);
    interval_right(right_node) = interval_right(i);
    
    arma::uvec left_idx, right_idx;
    if(interval_left(i) <= interval_right(left_node)-1){
      left_idx = arma::regspace<arma::uvec>(interval_left(i), interval_right(left_node)-1);
      count.col(left_node) = arma::sum(X.cols(left_idx), 1);
    }else{
      count.col(left_node) = arma::zeros<arma::vec>(X.n_rows);
    }
    
    
    if(interval_left(right_node) <= interval_right(i)-1){
      right_idx = arma::regspace<arma::uvec>(interval_left(right_node), interval_right(i)-1);
      count.col(right_node) = arma::sum(X.cols(right_idx), 1);
    }else{
      count.col(right_node) = arma::zeros<arma::vec>(X.n_rows);
    }
    
    
    arma::vec prob_i = count.col(left_node)/ count.col(i);
    prob_i.clamp(0.001, 0.999);
    prob_i.replace(arma::datum::nan, 0.5);
    
    empirical_phi.col(i) = arma::log(prob_i/(1-prob_i));
    
  }
  
  
  
  // see what this tree looks like
  Rcpp::List output;
  output = Rcpp::List::create(Rcpp::Named("count") = count,
                              Rcpp::Named("empirical_phi") = empirical_phi);
  
  return output;
}; 


