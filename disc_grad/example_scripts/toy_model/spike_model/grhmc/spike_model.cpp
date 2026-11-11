// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <random>
#include <vector>
#include <../../../../implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function.h>

struct spike_model { // plays the same role as model list in R code
  int n_parameters;
  double delta;
  double eta;
  
  Eigen::ArrayXd grad_indices;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd B;
  
  spike_model(int n_dim, double delta_value, double eta_value) // initialize the parameters
    : n_parameters(n_dim), delta(delta_value), eta(eta_value),
      grad_indices(Eigen::ArrayXd::Zero(n_dim)),
      A(3 * n_dim, n_dim)
  {
    // A = kronecker(diag(n), matrix(1, 3, 1))
    std::vector<Eigen::Triplet<double>> triplets; // Make the A matrix
    triplets.reserve(n_parameters * 3);
    for (int j = 0; j < n_parameters; ++j)
      for (int r = 0; r < 3; ++r)
        triplets.emplace_back(j * 3 + r, j, 1.0);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    
    // B = rep(c(-delta, 0, delta), n)
    Eigen::Vector3d pattern(-delta, 0.0, delta); // Make the B vector to give Az + B
    B = pattern.replicate(n_parameters, 1);
  }
  
  void grad_jump_fun(const Eigen::VectorXd& q) {
    if (q.size() != n_parameters) {
      throw std::invalid_argument("q does not match with n_parameters");
    }; // ensure that the size of grad_indices is the same as the number of parameters
    const auto qa = q.array();
    grad_indices = 1.0
    + (qa > -delta).cast<double>()
      + (qa > 0.0).cast<double>()
      + (qa > delta).cast<double>();
  }
  
  Eigen::VectorXd log_target_grad(const Eigen::VectorXd& q) {
    return -q.array()
    + (grad_indices == 2.0).cast<double>() * eta * (1.0 / delta)
    + (grad_indices == 3.0).cast<double>() * eta * (-1.0 / delta);
  }
};

//[[Rcpp::export]]
Rcpp::List grhmc_spike_model_cpp( // make an inner wrapper function that is input for the actual standard R wrapper function, but can run directly if correct list of settings is given
    Rcpp::List settings,
    double delta,
    double eta
) {
  
  spike_model spike_model_list(settings["n_parameters"], delta, eta);
  
  return grhmc_disc_grad::run_sampler(
    spike_model_list,
    settings
  );
  
}
