#include "PolyaGamma.h" // to sample polya-gamma random variable
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "../inst/include/CorTree_types.h"
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
  cortree::validate_counts(count_data);
  if (!tree.containsElementNamed("edge") || !tree.containsElementNamed("Nnode") ||
      !tree.containsElementNamed("tip.label")) {
    Rcpp::stop("tree must contain edge, Nnode, and tip.label.");
  }
  // Extract tree components using Armadillo types.
  arma::imat edge = Rcpp::as<arma::imat>(tree["edge"]); // two-column matrix: parent, child (1-indexed)
  int Nnode = Rcpp::as<int>(tree["Nnode"]);
  Rcpp::CharacterVector tip_label = tree["tip.label"];
  int n_tips = tip_label.size();
  int total_nodes = n_tips + Nnode;
  int n_samples = count_data.n_rows;
  arma::mat edge_numeric = Rcpp::as<arma::mat>(tree["edge"]);
  if (n_tips < 2 || Nnode != n_tips - 1 || count_data.n_cols != static_cast<arma::uword>(n_tips) ||
      edge.n_cols != 2 || edge.n_rows != static_cast<arma::uword>(total_nodes - 1) ||
      !edge_numeric.is_finite() || arma::any(arma::vectorise(edge_numeric) != arma::floor(arma::vectorise(edge_numeric))) ||
      arma::any(arma::vectorise(edge) < 1) || arma::any(arma::vectorise(edge) > total_nodes)) {
    Rcpp::stop("Require a rooted binary phylogeny with one count column per tip and valid node IDs.");
  }
  arma::ivec indegree(total_nodes, arma::fill::zeros);
  arma::ivec outdegree(total_nodes, arma::fill::zeros);
  for (arma::uword i = 0; i < edge.n_rows; ++i) {
    ++outdegree(edge(i,0) - 1);
    ++indegree(edge(i,1) - 1);
  }
  if (arma::accu(indegree == 0) != 1 || arma::any(indegree > 1) ||
      arma::any(outdegree.head(n_tips) != 0) || arma::any(outdegree.tail(Nnode) != 2)) {
    Rcpp::stop("tree must be rooted and binary, with each non-root node having one parent.");
  }
  arma::uvec roots = arma::find(indegree == 0);
  std::queue<int> pending;
  pending.push(static_cast<int>(roots(0)) + 1);
  int visited = 0;
  while (!pending.empty()) {
    int node = pending.front(); pending.pop(); ++visited;
    for (arma::uword i = 0; i < edge.n_rows; ++i)
      if (edge(i,0) == node) pending.push(edge(i,1));
  }
  if (visited != total_nodes) Rcpp::stop("tree must be connected and acyclic.");
  
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
    arma::mat parent_count;
    arma::mat kappa; // n(A) - n(A_p)/2, n by total_parents matrix

    arma::vec ind_layer_idx;
    arma::vec cor_layer_idx;
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
    arma::mat pi_full_sample;
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
  Rcpp::List test_output;
  bool all_ind = false;
  bool save_phi_trace = false;
  bool save_sigma_inv_trace = false;
  bool save_cluster_cor_trace = false;
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

  void initialize_mcmc_sample(){
    paras_sample.mu_sample = arma::zeros<arma::cube>(tree.total_parents, dat.n_clus, gibbs_control.mcmc_sample);
    if(save_phi_trace){
      paras_sample.phi_sample = arma::zeros<arma::cube>(dat.n, tree.total_parents, gibbs_control.mcmc_sample);
    }
    if(!all_ind && save_sigma_inv_trace){
      paras_sample.Sigma_inv_sample = Rcpp::List(gibbs_control.mcmc_sample);
    }
    paras_sample.pi_sample = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.mcmc_sample);
    paras_sample.pi_full_sample = arma::zeros<arma::mat>(dat.n_clus, gibbs_control.total_iter);
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
    tree.parent_count.set_size(dat.n, tree.total_parents);
    tree.kappa = arma::zeros<arma::mat>(dat.n, tree.total_parents);
    arma::uvec all_children = tree.edge.col(1);
    for(arma::uword j=0; j<tree.total_parents; j++){
      arma::uvec col_parent = find(tree.nodes == tree.parent_nodes(j));
      arma::uvec children = all_children(find(tree.edge.col(0) == tree.parent_nodes(j)));
      arma::uvec col_1st_child = find(tree.nodes == children(0));
      tree.parent_count.col(j) = tree.count.col(col_parent(0));
      tree.kappa.col(j) = tree.count.col(col_1st_child(0)) - tree.count.col(col_parent(0))/2;
    }

    // change this to be phylo-tree's layer
    tree.ind_layer_idx = 1.0 + parent_depth.elem(tree.idx_ind);
    tree.cor_layer_idx = 1.0 + parent_depth.elem(tree.idx_cor);
    

  }; 

  void set_test_output(){
    test_output = Rcpp::List::create(Rcpp::Named("count") = tree.count,
                                     Rcpp::Named("depth") = tree.m,
                                     Rcpp::Named("parent_count") = tree.parent_count,
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
    if (!all_ind && ghs_scale_hierarchy) {
      precision_scales.ones(dat.n_clus);
      precision_scale_rate = ghs_scale_shape > 1.0 ? ghs_scale_shape - 1.0 : ghs_scale_shape;
    }
    paras.mu = arma::zeros<arma::mat>(tree.total_parents, dat.n_clus);
    paras.Z_cor_status = arma::ones<arma::uvec>(dat.n_clus);
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
      arma::rowvec b_i = tree.parent_count.row(i);
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

    for (arma::uword i = 0; i < static_cast<arma::uword>(p); ++i) {
      arma::uvec ind = arma::regspace<arma::uvec>(0, p - 1);
      ind.shed_row(i);
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
    arma::uvec counts(dat.n_clus, arma::fill::zeros);
    for (arma::uword i = 0; i < paras.Z.n_elem; ++i) ++counts(paras.Z(i));
    paras.pi = cortree::stick_weights(counts, paras.alpha);
  };

  void update_loglike(){
    const arma::mat parent_count = tree.parent_count;
    paras.loglike = 0.0;
    for (arma::uword i = 0; i < parent_count.n_rows; ++i) {
      for (arma::uword j = 0; j < parent_count.n_cols; ++j) {
        double n = parent_count(i,j);
        paras.loglike += cortree::binomial_log_kernel(n, tree.kappa(i,j) + n / 2.0, paras.phi(i,j));
      }
    }

    if(save_cluster_cor_trace){
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
    paras_sample.pi_full_sample.col(iter) = paras.pi;
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
      if(save_cluster_cor_trace){
        paras_sample.cluster_cor.slice(idx) = paras.cluster_cor;
      }
      paras_sample.tausq_GHS.col(idx) = paras.tau_sq_GHS;
      if (!all_ind && ghs_scale_hierarchy) {
        paras_sample.scale_t.col(idx) = precision_scales;
        paras_sample.scale_b(idx) = precision_scale_rate;
      }
    }
    paras_sample.loglik(iter) = paras.loglike;
  };

  Rcpp::List get_gibbs_sample(){
    SEXP phi_output = save_phi_trace ? Rcpp::wrap(paras_sample.phi_sample) : R_NilValue;
    SEXP sigma_inv_output = save_sigma_inv_trace ? Rcpp::wrap(paras_sample.Sigma_inv_sample) : R_NilValue;
    SEXP cluster_cor_output = save_cluster_cor_trace ? Rcpp::wrap(paras_sample.cluster_cor) : R_NilValue;
    Rcpp::List output = Rcpp::List::create(Rcpp::Named("mu") = paras_sample.mu_sample,
                              Rcpp::Named("phi") = phi_output,
                              Rcpp::Named("Sigma_inv") = sigma_inv_output,
                              Rcpp::Named("cluster_cor") = cluster_cor_output,
                              Rcpp::Named("tausq_GHS") = paras_sample.tausq_GHS,
                              Rcpp::Named("pi_full") = paras_sample.pi_full_sample,
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
// [[Rcpp::export(signature = {count_data, tree, n_clus, cutoff_layer, total_iter, burnin, warm_start = 0L, init_Z, c_sigma2_vec = 1.0, sigma_mu2 = 1.0, all_ind = FALSE, cov_interval = 1L, save_phi_trace = FALSE, save_sigma_inv_trace = FALSE, save_cluster_cor_trace = FALSE, ghs_diag_rate = 0.0, ghs_diag_upper = Inf, ghs_jmlr_lambda = 0.0, ghs_det_df = 0.0, ghs_scale_hierarchy = FALSE, ghs_scale_shape = 3.0, ghs_scale_rate_shape = 2.0, ghs_scale_rate_rate = 1.0})]]
Rcpp::List PhyloTree_sampler(arma::mat count_data, Rcpp::List tree,
                      int n_clus, int cutoff_layer, 
                      int total_iter, int burnin, int warm_start=0,
                      arma::uvec init_Z = arma::zeros<arma::uvec>(1),
                      double c_sigma2_vec = 1.0, 
                      double sigma_mu2=1.0,
                      bool all_ind = false,
                      int cov_interval = 1,
                      bool save_phi_trace = false,
                      bool save_sigma_inv_trace = false,
                      bool save_cluster_cor_trace = false,
                      double ghs_diag_rate = 0.0,
                      double ghs_diag_upper = R_PosInf,
                      double ghs_jmlr_lambda = 0.0,
                      double ghs_det_df = 0.0,
                      bool ghs_scale_hierarchy = false,
                      double ghs_scale_shape = 3.0,
                      double ghs_scale_rate_shape = 2.0,
                      double ghs_scale_rate_rate = 1.0){
  Rcout<<"begin PhyloTree_sampler"<<std::endl;
  arma::wall_clock timer;
  timer.tic();
  cortree::validate_sampler(count_data, n_clus, total_iter, burnin, warm_start,
                            cov_interval, c_sigma2_vec, sigma_mu2, init_Z);
  cortree::validate_ghs_diag_rate(ghs_diag_rate);
  cortree::validate_ghs_jmlr_lambda(ghs_jmlr_lambda, ghs_diag_rate, ghs_diag_upper);
  cortree::validate_ghs_det_df(ghs_det_df, ghs_diag_upper, ghs_diag_rate, ghs_jmlr_lambda);
  cortree::validate_ghs_scale_hierarchy(ghs_scale_hierarchy, ghs_scale_shape,
    ghs_scale_rate_shape, ghs_scale_rate_rate, ghs_diag_rate, ghs_diag_upper,
    ghs_jmlr_lambda, warm_start);
  cortree::validate_ghs_diag_upper(ghs_diag_upper, ghs_diag_rate, warm_start, all_ind);
  if (cutoff_layer < 0) Rcpp::stop("cutoff_layer must be nonnegative.");
  PhyloTree model;

  if(init_Z.n_elem == 1 && count_data.n_rows > 1){
    init_Z = arma::randi<arma::uvec>(count_data.n_rows, arma::distr_param(0, n_clus-1));
  }
  
  model.all_ind = all_ind;
  model.save_phi_trace = save_phi_trace;
  model.save_sigma_inv_trace = save_sigma_inv_trace;
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
