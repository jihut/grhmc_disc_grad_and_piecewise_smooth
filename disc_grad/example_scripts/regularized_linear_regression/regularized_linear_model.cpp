// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <random>
#include <vector>
#include <../../implementation_scripts/general_scripts/cpp/grhmc_discontinuous_gradient_with_relevant_indices_function.h>

struct linear_model { // plays the same role as model list in R code
  int n_parameters;
  double mu;
  double rho;
  double gamma_mu;
  double gamma_sigma;
  Eigen::MatrixXd x;
  Eigen::VectorXd y;
  int length_beta;
  int n_rows;
  
  Eigen::ArrayXd grad_indices;
  Eigen::SparseMatrix<double> A;
  Eigen::VectorXd B;
  
  linear_model(
    int n_dim, 
    double mu_value, 
    double rho_value, 
    double gamma_mu_value, 
    double gamma_sigma_value,
    Eigen::MatrixXd x_matrix,
    Eigen::VectorXd y_vector
  ) // initialize the parameters
    : n_parameters(n_dim), mu(mu_value), rho(rho_value), gamma_mu(gamma_mu_value), gamma_sigma(gamma_sigma_value),
      x(x_matrix),
      y(y_vector),
      grad_indices(Eigen::ArrayXd::Zero(n_dim)),
      A(n_dim - 1, n_dim)
  {
    // A = kronecker(diag(n), matrix(1, 3, 1))
    for (int i = 0; i < n_dim - 1; i++) {
      A.insert(i, i) = 1;
    }
    
    B = Eigen::VectorXd::Zero(n_dim - 1);
    
    length_beta = (n_dim - 1) / 2;
    n_rows = x_matrix.rows();
  }
  
  void grad_jump_fun(const Eigen::VectorXd& q) {
    if (q.size() != n_parameters) {
      throw std::invalid_argument("q does not match with n_parameters");
    }; // ensure that the size of grad_indices is the same as the number of parameters
    for (int i = 0; i < n_parameters - 1; i++) {
      grad_indices(i) = (q(i) > 0) ? 1.0 : 0.0;
    }
  }
  
  Eigen::VectorXd log_target_grad(const Eigen::VectorXd& q) {
    Eigen::VectorXd beta = (q.array().segment(0, length_beta) * grad_indices.array().segment(0, length_beta)) - (q.array().segment(length_beta, length_beta) * grad_indices.array().segment(length_beta, length_beta));
    Eigen::VectorXd y_minus_linear_predictor = y - x * beta; 
    double exp_gamma = exp(q[n_parameters - 1]);
    Eigen::VectorXd ratio_grad_indices_exp_gamma = grad_indices / exp_gamma; 
    Eigen::VectorXd log_target_grad_value_vec(n_parameters);
    Eigen::VectorXd prior_beta_minus_plus_term = -(q.array() - mu) / (std::pow(rho, 2));
    Eigen::VectorXd y_minus_linear_predictor_times_scale_x = y_minus_linear_predictor.transpose() * x; 
    // If grad_indices = 0 for i = 1:length_beta --> beta_i+ < 0 so only the prior part in this region, no term related to observations and covariates
    log_target_grad_value_vec.segment(0, length_beta) = 
      prior_beta_minus_plus_term.array().segment(0, length_beta) + 
      y_minus_linear_predictor_times_scale_x.array() * ratio_grad_indices_exp_gamma.array().segment(0, length_beta);
    log_target_grad_value_vec.segment(length_beta, length_beta) = 
      prior_beta_minus_plus_term.array().segment(length_beta, length_beta) - 
      y_minus_linear_predictor_times_scale_x.array() * ratio_grad_indices_exp_gamma.array().segment(length_beta, length_beta);
    log_target_grad_value_vec(n_parameters - 1) = 
      -(n_rows / 2.0) - (q[n_parameters - 1] - gamma_mu) / (gamma_sigma * gamma_sigma) + (y_minus_linear_predictor.array() * y_minus_linear_predictor.array()).sum() / (2.0 * exp_gamma); 
      
      return log_target_grad_value_vec; 
  }
};

//[[Rcpp::export]]
Rcpp::List grhmc_linear_model_cpp( // make an inner wrapper function that is input for the actual standard R wrapper function, but can run directly if correct list of settings is given
    Rcpp::List settings,
    double mu, 
    double rho, 
    double gamma_mu, 
    double gamma_sigma,
    Eigen::MatrixXd x_matrix,
    Eigen::VectorXd y_vector
) {
  
  linear_model linear_model_list(
      settings["n_parameters"], 
              mu, 
              rho,
              gamma_mu,
              gamma_sigma,
              x_matrix,
              y_vector);
  
  return grhmc_disc_grad::run_sampler(
    linear_model_list,
    settings
  );
  
}
