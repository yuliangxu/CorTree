#include "PolyaGamma.h" // to sample polya-gamma random variable
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <queue>
#include <algorithm>
#include <functional>


//' Aggregate counts along a phylogenetic tree and output node depths
 //'
 //' This function takes a count matrix (with rows as samples and columns corresponding to tip taxa)
 //' along with a phylo tree (an R list with elements "edge", "Nnode", "tip.label") and returns
 //' a list with two items. The first item, "aggregated", is an aggregated count matrix where each row is a sample,
 //' the first column holds the counts at the root (i.e. total counts), and subsequent columns hold counts for nodes deeper in the tree.
 //' The second item, "depth", is an integer vector giving the depth (distance from the root)
 //' for each node corresponding to the columns of the aggregated matrix.
 //'
 //' @param count_data An n x m matrix of counts (n samples, m tips)
 //' @param tree A list representing a phylo tree (must contain "edge", "Nnode", "tip.label")
 //' @return A list with two elements: "aggregated" (the aggregated count matrix) and "depth" (the depth for each node).
 //' @export
 // [[Rcpp::export]]
 Rcpp::List aggregate_tree_counts(arma::mat count_data, Rcpp::List tree) {
  // Extract tree components using Armadillo types.
  arma::imat edge = Rcpp::as<arma::imat>(tree["edge"]); // two-column matrix: parent, child (1-indexed)
  int Nnode = Rcpp::as<int>(tree["Nnode"]);
  Rcpp::CharacterVector tip_label = tree["tip.label"];
  int n_tips = tip_label.size();
  int total_nodes = n_tips + Nnode;
  int n_samples = count_data.n_rows;
  
  // Create an aggregated count matrix.
  // For tip nodes (columns 0 to n_tips-1) we use the original count_data.
  arma::mat agg(n_samples, total_nodes, arma::fill::zeros);
  agg.cols(0, n_tips - 1) = count_data;
  
  // Build a children list for each node using an Armadillo field.
  // First, count how many children each node has.
  arma::Col<int> children_count(total_nodes, arma::fill::zeros);
  for (arma::uword i = 0; i < edge.n_rows; i++) {
    int parent = edge(i, 0);
    children_count(parent - 1)++;  // adjust from 1-indexed to 0-indexed
  }
  
  // Allocate the children field: each entry is a vector of children for that node.
  arma::field<arma::Col<int>> children(total_nodes);
  for (int i = 0; i < total_nodes; i++) {
    children(i) = arma::Col<int>(children_count(i));
  }
  
  // Fill the children field and mark nodes that are children.
  arma::Col<int> child_index(total_nodes, arma::fill::zeros);
  arma::Col<int> isChild(total_nodes, arma::fill::zeros); // 0 means not a child, 1 means child
  
  for (arma::uword i = 0; i < edge.n_rows; i++) {
    int parent = edge(i, 0);
    int child = edge(i, 1);
    children(parent - 1)(child_index(parent - 1)) = child;
    child_index(parent - 1)++;
    isChild(child - 1) = 1;
  }
  
  // Identify the root as the node that is never marked as a child.
  int root = -1;
  for (int i = 0; i < total_nodes; i++) {
    if (isChild(i) == 0) {
      root = i + 1;  // convert back to 1-indexed
      break;
    }
  }
  if (root == -1) {
    Rcpp::stop("Could not find root in the tree.");
  }
  
  // Recursive aggregation function: for internal nodes, add counts from each child.
  std::function<void(int)> recurse = [&](int node) {
    if (node <= n_tips) return; // tip node: counts are already set
    arma::Col<int>& child_list = children(node - 1);
    for (arma::uword j = 0; j < child_list.n_elem; j++) {
      int child = child_list(j);
      recurse(child);
      agg.col(node - 1) += agg.col(child - 1);
    }
  };
  
  // Start aggregation from the root.
  recurse(root);
  
  // Compute the depth (distance from the root) for every node using an Armadillo integer vector.
  arma::ivec depth(total_nodes, arma::fill::value(-1));
  depth(root - 1) = 0;
  std::queue<int> q;
  q.push(root);
  while (!q.empty()) {
    int cur = q.front();
    q.pop();
    arma::Col<int>& child_list = children(cur - 1);
    for (arma::uword j = 0; j < child_list.n_elem; j++) {
      int child = child_list(j);
      depth(child - 1) = depth(cur - 1) + 1;
      q.push(child);
    }
  }
  
  // Create a vector of node indices (1-indexed) and sort them by their depth.
  std::vector<int> nodes(total_nodes);
  for (int i = 0; i < total_nodes; i++) {
    nodes[i] = i + 1;
  }
  std::sort(nodes.begin(), nodes.end(), [&](int a, int b) {
    if (depth(a - 1) == depth(b - 1))
      return a < b;
    return depth(a - 1) < depth(b - 1);
  });
  
  // Build the output matrix with columns arranged by increasing depth.
  arma::mat result(n_samples, total_nodes);
  arma::ivec sorted_depth(total_nodes);
  for (int j = 0; j < total_nodes; j++) {
    int node = nodes[j];
    result.col(j) = agg.col(node - 1);
    sorted_depth(j) = depth(node - 1);
  }
  
  // Generate the parent_nodes vector:
  // Iterate over the sorted nodes and keep only those with at least one child.
  std::vector<int> parent_nodes;
  for (int j = 0; j < total_nodes; j++) {
    int node = nodes[j];
    if (children_count(node - 1) > 0) {  // node has children → it is a parent node
      parent_nodes.push_back(node);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("aggregated") = result,
                            Rcpp::Named("nodes") = nodes,
                            Rcpp::Named("parent_nodes") = parent_nodes,
                            Rcpp::Named("depth") = sorted_depth);
}

// [[Rcpp::export]]
arma::uvec isIn(const arma::uvec& a, const arma::uvec& b) {
  arma::uvec result(a.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < a.n_elem; i++) {
    // Check membership: if a(i) is in b, then set result(i) to 1.
    if (arma::any(b == a(i))) {
      result(i) = 1;
    }
  }
  return result;
}

// [[Rcpp::export]]
arma::uvec complementarySet(const arma::uvec& a, const arma::uvec& b) {
  // We use a std::vector to collect the complementary elements
  std::vector<arma::uword> comp;
  
  // Loop over each element in b.
  for (arma::uword i = 0; i < b.n_elem; i++) {
    // If b[i] is not found in a, then include it.
    if (!arma::any(a == b[i])) {
      comp.push_back(b[i]);
    }
  }
  
  // Convert the std::vector to an arma::uvec and return.
  return arma::uvec(comp);
}

class PhyloTree{
private: 
  struct Data{
    arma::mat X; // n by (n_species) matrix, count table, sparse matrix
    
    int n_clus; // number of latent clusters 
    int n; // number of observations n
  } dat;
  
  struct Tree{
    int m; // tree depth. layer m=0 is the sample space Omega
    
    int cutoff_layer; // cutoff layer for correlated nodes, 0,...,m
    int L; // cutoff for the number of correlated nodes in a tree, any number in {2^{cutoff_layer+1}-1, i=0,1,...,m}
    arma::uvec parent_nodes; // parent nodes, sorted by depth
    arma::uvec nodes; // all nodes, sorted by depth
    arma::umat edge; // edge matrix, parent, child

    arma::uvec idx_cor; // length of L
    arma::uvec idx_ind; // length of (total_parents-L)
    // for mu, fixed mu for the last layer
    arma::uvec idx_ind_mu; // length of (total_parents-L)
    arma::uvec idx_fixed_mu; // index of independent nodes in the last layer
    
    int total_nodes; // total number of nodes, = 2^{m+1}-1
    int total_parents; // total number of parents, = 2^m-1, also the total number of phi
    arma::mat count; // X_i(A_eps), n by total_nodes matrix, count table 
    arma::mat kappa; // n(A) - n(A_p)/2, n by total_parents matrix

    arma::vec ind_layer_idx;
  } tree; // fixed tree structure
  
  struct TreeLatentParameters{
    arma::mat mu; //  total_parents x n_clus, node probability
    arma::cube Sigma_inv; // L by L x n_clus precision matrix
    arma::mat sigma2_vec; // length of (total_parents-L) by n_clus matrix, diagonal elements of independent nodes. fixed
    arma::mat omega; // n by total_parents matrix, polya gamma latent variable
    
    // Dirichlet process  mixure parameters
    arma::uvec Z; // n by 1, membership indicator {1,2,...,n_clus}
    arma::uvec Z_cor_status; // n by 1, membership indicator {1,2,...,n_clus}
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
  bool save_phi_trace = false;
  bool save_sigma_inv_trace = false;
  int warm_start = 0;
  int cov_interval = 1;

  
  
  
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
    if(save_phi_trace){
      paras_sample.phi_sample = arma::zeros<arma::cube>(dat.n, tree.total_parents, gibbs_control.mcmc_sample);
    }
    if(!all_ind && save_sigma_inv_trace){
      paras_sample.Sigma_inv_sample = Rcpp::List(gibbs_control.mcmc_sample);
    }
    paras_sample.pi_sample = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.mcmc_sample);
    paras_sample.loglik = arma::zeros<arma::vec>(gibbs_control.total_iter);
    paras_sample.Z_sample = arma::zeros<arma::umat>(dat.n, gibbs_control.mcmc_sample);
    
    paras_sample.cluster_cor = arma::zeros<arma::ucube>(dat.n, dat.n, gibbs_control.mcmc_sample);
    paras_sample.tausq_GHS = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.mcmc_sample);
  };
  
  
  void vectorize_tree(arma::mat count_data, Rcpp::List phylo_tree, arma::uword cutoff_layer){
    Rcpp::List aggregate_tree = aggregate_tree_counts(count_data, phylo_tree);
    arma::vec node_depth = aggregate_tree["depth"];
    
    int Nnode = Rcpp::as<int>(phylo_tree["Nnode"]);
    arma::uvec parent_nodes = aggregate_tree["parent_nodes"];
    arma::umat edge = phylo_tree["edge"];
    tree.edge = edge;
    
    // set tree structure
    tree.m = max(node_depth);
    tree.total_nodes = node_depth.size();
    tree.total_parents = Nnode;
    tree.parent_nodes = parent_nodes;
    arma::uvec all_nodes = aggregate_tree["nodes"];
    tree.nodes = all_nodes;
    
    if(cutoff_layer > tree.m){
      Rcpp::stop("cutoff_layer should be less than the tree depth");
    }else{
      tree.cutoff_layer = cutoff_layer;
    }
    
    arma::vec parent_depth(tree.total_parents, arma::fill::zeros);
    for(arma::uword j = 0; j < tree.total_parents; j++){
      arma::uvec node_pos = arma::find(tree.nodes == tree.parent_nodes(j), 1, "first");
      if(node_pos.n_elem == 0){
        Rcpp::stop("Parent node not found in sorted node list.");
      }
      parent_depth(j) = node_depth(node_pos(0));
    }

    if(all_ind){
      tree.L = 0;
      tree.idx_ind = arma::regspace<arma::uvec>(0, tree.total_parents-1);
    }else{
      // correlated parent indices use parent-parameter indexing (0..total_parents-1)
      tree.idx_cor = arma::find(parent_depth <= static_cast<double>(cutoff_layer));
      tree.L = tree.idx_cor.n_elem;
      tree.idx_ind = complementarySet(tree.idx_cor, arma::regspace<arma::uvec>(0, tree.total_parents-1));
    }
    

    // ser precision for Sigma_inv
    hyper.err_precision = (tree.L > 0) ? pow(10.0, -11.0/sqrt(static_cast<double>(tree.L))) : 1e-11;
    Rcpp::Rcout<<"err_precision="<<hyper.err_precision<<std::endl;
    
    arma::mat aggregated = aggregate_tree["aggregated"];
    tree.count = aggregated;

    // initialize kappa as n(A) - n(A_p)/2, columns are ordered the same way as the nodes
    tree.kappa = arma::zeros<arma::mat>(dat.n, tree.total_parents);
    arma::uvec all_children = tree.edge.col(1);
    for(arma::uword j=0; j<tree.total_parents; j++){
      arma::uvec col_parent = find(tree.nodes == tree.parent_nodes(j));
      arma::uvec children = all_children(find(tree.edge.col(0) == tree.parent_nodes(j)));
      arma::uvec col_1st_child = find(tree.nodes == children(0));
      tree.kappa.col(j) = tree.count.col(col_1st_child(0)) - tree.count.col(col_parent(0))/2;
    }

    // change this to be phylo-tree's layer
    tree.ind_layer_idx = parent_depth.elem(tree.idx_ind);
    

  }; 

  void set_test_output(){
    test_output = Rcpp::List::create(Rcpp::Named("count") = tree.count,
                                     Rcpp::Named("depth") = tree.m,
                                     Rcpp::Named("tree.idx_ind") = tree.idx_ind,
                                     Rcpp::Named("tree.idx_cor") = tree.idx_cor,
                                     Rcpp::Named("ind_layer_idx") = tree.ind_layer_idx,
                                     Rcpp::Named("sigma_vec") = paras.sigma2_vec,
                                     Rcpp::Named("kappa") = tree.kappa);
  }
  

  // arma::uword node_to_column(arma::uword node){
  //   arma::uvec indices = arma::find(tree.nodes == node);
  //   return indices(0);
  // };
  
  // arma::uword column_to_node(arma::uword column){
  //   return tree.nodes(column);
  // };

  void initialize_parameters(arma::uvec init_Z){
    paras.mu = arma::zeros<arma::mat>(tree.total_parents, dat.n_clus);
    paras.Z_cor_status = arma::ones<arma::uvec>(dat.n);
    // Assuming tree.L is the dimension of the identity matrix and tree.total_slices is the number of slices

    if(all_ind){
      paras.sigma2_vec = arma::ones<arma::mat>(tree.total_parents, dat.n_clus);
    }else{
      paras.Sigma_inv = arma::cube(tree.L, tree.L, dat.n_clus, arma::fill::zeros);
      for (size_t i = 0; i < dat.n_clus; i++) {
          paras.Sigma_inv.slice(i) = arma::eye<arma::mat>(tree.L, tree.L);
      }
      paras.sigma2_vec = arma::ones<arma::mat>(tree.total_parents - tree.L, dat.n_clus);
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
      // Old wrong code: indexed phi by cluster label instead of sample index.
      // arma::rowvec c_i = paras.phi.row(paras.Z(i));
      arma::rowvec c_i = paras.phi.row(i);
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
    // Rcout<<"begin of GHS_oneSample"<<std::endl;
    // Rcout<<"is S pd?"<<S_pd_check<<std::endl;
    // Rcout<<"is Omega pd?"<<pd_check<<std::endl;

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
            
            // update sigma_vec
            arma::vec layer_idx = floor( log2(2+arma::regspace<arma::vec>(0, tree.L - 1)) );
            arma::vec a_vec = layer_idx* hyper.c_sigma2_vec + idx_k.n_elem/2;
            arma::mat phi_res = paras.phi(idx_k,tree.idx_cor);
            phi_res.each_row() -= paras.mu(tree.idx_cor, arma::uvec{k}).t();
            arma::vec b_vec = layer_idx + arma::trans(arma::sum(phi_res%phi_res, 0)/2);  
            
            arma::vec gamma1(a_vec.n_elem);
            arma::vec gamma2(b_vec.n_elem);
            for (arma::uword i = 0; i < a_vec.n_elem; ++i) {
              gamma1(i) = arma::randg(arma::distr_param(a_vec(i), 1.0));  // generate gamma random values
              gamma2(i) = arma::randg(arma::distr_param(b_vec(i), 1.0));
            }
            arma::vec sigma2_k_vec = (gamma1 + gamma2)/gamma1;
            
            
            paras.Sigma_inv.slice(k) = arma::diagmat(sigma2_k_vec);
            
            
          }else if(arma::det(paras.Sigma_inv.slice(k)) < 1e200){
            
            arma::mat phi_centered = paras.phi(idx_k, tree.idx_cor);
            phi_centered.each_row() -= paras.mu(tree.idx_cor, arma::uvec{k}).t();
            // arma::rowvec phi_mean = arma::mean(phi_centered, 0);
            // phi_centered.each_row() -= phi_mean;
            arma::mat phi_centered_Cov = phi_centered.t() * phi_centered;
            GHS_oneSample(phi_centered_Cov, idx_k.n_elem, ghs_list[k]);
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
        arma::vec a_vec = tree.ind_layer_idx * hyper.c_sigma2_vec + idx_k.n_elem/2;
        arma::mat phi_res = paras.phi(idx_k,tree.idx_ind);
        phi_res.each_row() -= paras.mu(tree.idx_ind, arma::uvec{k}).t();
        arma::vec b_vec = tree.ind_layer_idx + arma::trans(arma::sum(phi_res%phi_res, 0)/2);  
        
        arma::vec gamma1(a_vec.n_elem);
        arma::vec gamma2(b_vec.n_elem);
        for (arma::uword i = 0; i < a_vec.n_elem; ++i) {
          gamma1(i) = arma::randg(arma::distr_param(a_vec(i), 1.0));  // generate gamma random values
          gamma2(i) = arma::randg(arma::distr_param(b_vec(i), 1.0));
        }
        paras.sigma2_vec.col(k) = (gamma1 + gamma2)/gamma1;
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
          // Old wrong code: missing cluster-dependent Gaussian normalization
          // term (+0.5 * log|Sigma_inv_k|), which can bias cluster assignment.
          // loglik_cor_i(k) = -0.5*arma::dot(diff_sub, paras.Sigma_inv.slice(k) * diff_sub);
          double log_det_val = 0.0;
          double sign_det = 0.0;
          arma::log_det(log_det_val, sign_det, paras.Sigma_inv.slice(k));
          if(sign_det <= 0.0 || !arma::is_finite(log_det_val)){
            loglik_cor_i(k) = -arma::datum::inf;
          }else{
            loglik_cor_i(k) = -0.5 * arma::dot(diff_sub, paras.Sigma_inv.slice(k) * diff_sub) + 0.5 * log_det_val;
          }

        }

        arma::vec mu_ind = paras.mu(tree.idx_ind, arma::uvec{k});
        arma::vec diff_ind = phi_ind - mu_ind;
        // Old wrong code: missing diagonal-Gaussian normalization term
        // (-0.5 * sum(log(sigma2_k))), which can bias cluster assignment.
        // loglik_ind_i(k) = -0.5*arma::dot(diff_ind, diff_ind/paras.sigma2_vec.col(k));
        arma::vec sigma2_k = paras.sigma2_vec.col(k);
        if(arma::any(sigma2_k <= 0.0) || !sigma2_k.is_finite()){
          loglik_ind_i(k) = -arma::datum::inf;
        }else{
          double quad_ind = arma::dot(diff_ind, diff_ind / sigma2_k);
          double log_det_ind = arma::sum(arma::log(sigma2_k));
          loglik_ind_i(k) = -0.5 * (quad_ind + log_det_ind);
        }
      }
      
      arma::vec loglik_i = loglik_ind_i + loglik_cor_i;
      
      arma::vec log_prob = loglik_i + log(paras.pi);
      // Old wrong code: if all entries are non-finite (e.g., all -Inf), this
      // normalization creates NaN probabilities and invalid CDF sampling.
      // log_prob -= arma::max(log_prob);
      // log_prob = arma::exp(log_prob);
      // log_prob /= arma::sum(log_prob);
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
    if(iter >= gibbs_control.burnin){
      int idx = iter - gibbs_control.burnin;
      paras_sample.mu_sample.slice(idx) = paras.mu;
      if(!all_ind && save_sigma_inv_trace){
        paras_sample.Sigma_inv_sample[idx] = paras.Sigma_inv;
      }
      paras_sample.pi_sample.col(idx) = paras.pi;
      if(save_phi_trace){
        paras_sample.phi_sample.slice(idx) = paras.phi;
      }
      paras_sample.Z_sample.col(idx) = paras.Z;
      paras_sample.cluster_cor.slice(idx) = paras.cluster_cor;
      paras_sample.tausq_GHS.col(idx) = paras.tau_sq_GHS;
    }
    paras_sample.loglik(iter) = paras.loglike;
  };

  Rcpp::List get_gibbs_sample(){
    SEXP phi_output = save_phi_trace ? Rcpp::wrap(paras_sample.phi_sample) : R_NilValue;
    SEXP sigma_inv_output = save_sigma_inv_trace ? Rcpp::wrap(paras_sample.Sigma_inv_sample) : R_NilValue;
    return Rcpp::List::create(Rcpp::Named("mu") = paras_sample.mu_sample,
                              Rcpp::Named("phi") = phi_output,
                              Rcpp::Named("Sigma_inv") = sigma_inv_output,
                              Rcpp::Named("cluster_cor") = paras_sample.cluster_cor,
                              Rcpp::Named("tausq_GHS") = paras_sample.tausq_GHS,
                              Rcpp::Named("pi") = paras_sample.pi_sample,
                              Rcpp::Named("Z") = paras_sample.Z_sample,
                              Rcpp::Named("loglik") = paras_sample.loglik);
  };
  
  
};

// [[Rcpp::export]]
Rcpp::List PhyloTree_sampler(arma::mat count_data, Rcpp::List tree,
                      int n_clus, int cutoff_layer, 
                      int total_iter, int burnin, int warm_start=0,
                      arma::uvec init_Z = arma::zeros<arma::uvec>(1),
                      double c_sigma2_vec = 1.0, 
                      double sigma_mu2=1.0,
                      bool all_ind = false,
                      int cov_interval = 1,
                      bool save_phi_trace = false,
                      bool save_sigma_inv_trace = false){
  Rcout<<"begin PhyloTree_sampler"<<std::endl;
  arma::wall_clock timer;
  timer.tic();
  PhyloTree model;

  if(init_Z.n_elem == 1){
    init_Z = arma::randi<arma::uvec>(count_data.n_rows, arma::distr_param(0, n_clus-1));
  }
  
  model.all_ind = all_ind;
  model.save_phi_trace = save_phi_trace;
  model.save_sigma_inv_trace = save_sigma_inv_trace;
  model.warm_start = warm_start;
  model.cov_interval = cov_interval;
  model.load_data(count_data, n_clus);
  Rcpp::Rcout << "Data loaded" << std::endl;
  model.set_hyperparameter(c_sigma2_vec, sigma_mu2);
  Rcpp::Rcout << "Hyperparameter set" << std::endl;
  model.set_gibbs_control(total_iter, burnin);
  Rcpp::Rcout << "Gibbs control set" << std::endl;
  model.vectorize_tree(count_data, tree, cutoff_layer);
  Rcpp::Rcout << "Tree vectorized" << std::endl;
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
