#include <RcppEigen.h>
#include <random>
//[[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include <cmath>
#include <vector>
// #include <boost/math/tools/cubic_roots.hpp>
#include <unsupported/Eigen/Polynomials>
#include <utility>
#include <chrono>


namespace grhmc_disc_grad {

// A model list similar to R that should contain grad_jump_fun, the matrix A and B defining the linear boundary equations Aq + B,
// initializing grad_indices and containing the gradient. 
template<class model_list>
struct model_spec_fun{ // contains all the functions that are induced by the model given by the model list that will be used in running the sampler
  
  model_list& model_list_obj; // initialize a model list object
  int n_evals_ode = 0;
  double lambda;   
  int n_parameters;
  int n_root_eqs;
  double push_eps_tol;
  int n_boundary_crossings = 0;
  std::vector<std::vector<int>> non_zero_cols_of_A; // for each row of A, store the non-zero column indices
  Eigen::VectorXd m_elements;
  Eigen::VectorXd s_elements;
  model_spec_fun( // arguments to initialize a model_spec_fun object
    model_list& model_list_init,
    double lambda_init,
    int n_parameters_init,
    double push_eps_tol_init,
    Eigen::VectorXd m_elements_init,
    Eigen::VectorXd s_elements_init
  ) : model_list_obj(model_list_init), // pass the given model list object to initialize this struct containing all the relevant functions baed on the model list
  lambda(lambda_init), n_parameters(n_parameters_init),
  push_eps_tol(push_eps_tol_init),
  m_elements(m_elements_init),
  s_elements(s_elements_init)
  {
    n_root_eqs = 1 + (int)model_list_obj.A.rows();
    non_zero_cols_of_A.resize(model_list_obj.A.rows());
    for (int i = 0; i < model_list_obj.A.rows(); ++i) {
      for (int j = 0; j < model_list_obj.A.cols(); ++j) {
        if (model_list_obj.A.coeff(i, j) != 0.0) {
          non_zero_cols_of_A[i].push_back(j);  
        }
      }
    }
  }
  
  
  Eigen::VectorXd ode( // ODE itself
      const double& t, 
      const Eigen::VectorXd& state
  ) {
    ++n_evals_ode; 
    Eigen::VectorXd ret(2 * n_parameters + 3);
    ret(0) = 0; // number of momentum refresh event
    ret(1) = lambda; // integrate lambda to get Lambda
    ret(2) = 0.0; // u
    ret.segment(3, n_parameters) = state.segment(3 + n_parameters, n_parameters); // \dot qbar = pbar
    ret.segment(3 + n_parameters, n_parameters) = model_list_obj.log_target_grad(
      m_elements + state.segment(3, n_parameters).cwiseProduct(s_elements)
    ).cwiseProduct(s_elements);
    
    return ret;
    
  }
  
  Eigen::VectorXd lin_root_fun(
      const double& t,
      const Eigen::VectorXd& state
  ) {
    Eigen::VectorXd ret(n_root_eqs);
    ret(0) = state(1) - state(2); // momentum refresh
    ret.segment(1, n_root_eqs - 1) = model_list_obj.A * (state.segment(3, n_parameters).cwiseProduct(s_elements) + m_elements)
      + model_list_obj.B; // boundary crossing equations
    return ret;
  }
  
  Eigen::VectorXd region_lin_event_fun(
      const double& t,
      const Eigen::VectorXd& state
  ) {
    Eigen::VectorXd abs_yroot = lin_root_fun(t, state).array().abs();
    double min_abs_yroot = abs_yroot.minCoeff();
    int index_min_abs_yroot; 
    
    for (int i = 0; i < n_root_eqs; i++) {
      if (abs_yroot[i] == min_abs_yroot) {
        index_min_abs_yroot = i; 
      }
    }
    Eigen::VectorXd new_state = state; 
    
    if (index_min_abs_yroot == 0) { // momentum refresh event
      // std::cout << "momentum refresh at t: " << t << std::endl;
      new_state(0) = state(0) + 1; 
      new_state(1) = 0.0; 
      new_state(2) = R::exp_rand();
      // state.segment(3, n_parameters) = state.segment(3, n_parameters);
      for (int i = 0; i < n_parameters; ++i) {
        new_state(3 + n_parameters + i) = R::norm_rand();
      }
    } else { // root due to crossing boundaries of discontinuous gradient
      Eigen::VectorXd qbar = state.segment(3, n_parameters);
      Eigen::VectorXd pbar = state.segment(3 + n_parameters, n_parameters);
      // Eigen::VectorXd qp = q.array() + 0.0001 * p.array();
      ++n_boundary_crossings;
      Eigen::VectorXd qp = qbar;
      int num_of_non_zero_columns = non_zero_cols_of_A[index_min_abs_yroot - 1].size();
      for (int j = 0; j < num_of_non_zero_columns; j++) { // only shift the coordinates affected by the root
        qp(non_zero_cols_of_A[index_min_abs_yroot - 1][j]) = qbar(non_zero_cols_of_A[index_min_abs_yroot - 1][j]) + push_eps_tol * pbar(non_zero_cols_of_A[index_min_abs_yroot - 1][j]);
      }
      
      model_list_obj.grad_jump_fun(qp.cwiseProduct(s_elements) + m_elements); // update grad_indices
      new_state.segment(0, 3) = state.segment(0, 3);
      new_state.segment(3, n_parameters) = qbar; 
      new_state.segment(3 + n_parameters, n_parameters) = pbar; 
    }
    
    return new_state;
  }
  
};

struct root_info { // make a struct to store info about previous and current root
  int which_root = -1;
  double root_time = 2;
  int root_type = -1;
};

Eigen::VectorXd CubicPolyRealRoots(const double& a, const double& b, const double& c, const double& d, const double& precision_real_root_lin_root_finder) { // not used anymore
  
  if (std::fabs(a) < precision_real_root_lin_root_finder) { // if this is true --> the root function is at most a quadratic polynomial
    
    if (std::fabs(b) < precision_real_root_lin_root_finder) { // if this is true --> the root function is at most linear
      
      if (std::fabs(c) < precision_real_root_lin_root_finder) { // if this is true --> the root function is constant
        // std::cout << "Constant function!" << std::endl; 
        Eigen::VectorXd ret(0);
        return ret; 
        
      } else {
        
        // std::cout << "Linear polynomial!" << std::endl; 
        
        Eigen::VectorXd ret(1);
        ret(0) = -d / c; 
        
        return ret; 
        
      }
      
    } else {
      
      // std::cout << "Quadratic polynomial!" << std::endl; 
      
      // ABC formula 
      
      double quadratic_discriminant = std::pow(c, 2) - 4.0 * b * d;
      
      if (quadratic_discriminant >= 0.0) {
        
        Eigen::VectorXd ret(2);
        ret(0) = (-c + std::sqrt(quadratic_discriminant)) / (2.0 * b);
        ret(1) = (-c - std::sqrt(quadratic_discriminant)) / (2.0 * b);
        
        return ret; 
        
      } else {
        
        Eigen::VectorXd ret(0); 
        // ret(0) = std::numeric_limits<double>::max(); // dummy to have something to return when no real roots
        
        return ret; 
        
      }
      
    }
    
    
  } else {
    // std::cout << "Cubic polynomial!" << std::endl; 
    // Using Cardano's formula
    
    double Q = (3.0 * a * c - std::pow(b, 2)) / (9.0 * std::pow(a, 2));
    double R = (9.0 * a * b * c - 27.0 * std::pow(a, 2) * d - 2.0 * std::pow(b, 3)) / (54.0 * std::pow(a, 3));
    double QQQ = std::pow(Q, 3);
    double RR = std::pow(R, 2);
    double D = QQQ + RR;
    
    if (D < 0) { // If discriminant is less than zero --> three distinct real roots
      double theta = std::acos(R / (std::sqrt(-QQQ)));
      
      Eigen::VectorXd ret(3); 
      ret[0] = 2.0 * std::sqrt(-Q) * std::cos(theta / 3.0) - b / (3.0 * a);
      ret[1] = 2.0 * std::sqrt(-Q) * std::cos(theta / 3.0 + (2.0 * M_PI) / 3.0) - b / (3.0 * a);
      ret[2] = 2.0 * std::sqrt(-Q) * std::cos(theta / 3.0 + (4.0 * M_PI) / 3.0) - b / (3.0 * a);
      
      return ret;
      
    } else { // otherwise one real root (assuming that the discriminant is non-zero)
      
      double sqrtD = std::sqrt(D);
      double S = std::cbrt(R + sqrtD);
      double T = std::cbrt(R - sqrtD);
      
      Eigen::VectorXd ret(1);
      ret[0] = S + T  - b / (3.0 * a);
      
      return ret;
      
    }
  }
  
}

Eigen::PolynomialSolver<double, Eigen::Dynamic> solver; // rather use Polynomial solver of Eigen

Eigen::VectorXd determine_poly_degree( // ax^3 + bx^2 + cx + d
    double a,
    double b,
    double c,
    double d,
    double coeff_tol // the tolerance such that we define a real root if abs(Imaginary part) < 1e-13 similar to R function
) {
  
  if (std::fabs(a) < coeff_tol) { // second degree polynomial at most
    if (std::fabs(b) < coeff_tol) { // first degree polynomial at most
      if (std::fabs(c) < coeff_tol) { // constant function
        Eigen::VectorXd ret(0);
        return ret;
      } 
      Eigen::VectorXd ret(2);
      ret(0) = d;
      ret(1) = c;
      return ret;
    }
    Eigen::VectorXd ret(3);
    ret(0) = d;
    ret(1) = c;
    ret(2) = b;
    return ret;
  }
  Eigen::VectorXd ret(4);
  ret(0) = d;
  ret(1) = c;
  ret(2) = b;
  ret(3) = a;
  return ret;
  
}

// Construct Bernstein coefficients and check if the sign of all the coefficients matches
// If not, then one knows that there is no root between 0 and 1 and can drop the solver step
bool cubic_root_check( // ax^3 + bx^2 + cx + d
    double a,
    double b,
    double c,
    double d
) {
  
  double b0 = d;
  double b1 = d + c / 3;
  double b2 = d + 2 * c / 3 + b / 3;
  double b3 = a + b + c + d;
  
  return (b0 > 0.0 && b1 > 0.0 && b2 > 0.0 && b3 > 0.0) ||
    (b0 < 0.0 && b1 < 0.0 && b2 < 0.0 && b3 < 0.0);
  
}

template<class model_spec>
root_info deSolverLinRootFinder(
    model_spec& model_spec_obj,  
    const Eigen::VectorXd& state,
    const Eigen::VectorXd& new_state,
    const Eigen::VectorXd& k1,
    const Eigen::VectorXd& k4, 
    const double t,
    const double h,
    const root_info& last_root_obj,
    const double last_root_offset,
    const double precision_real_root_finder,
    const double coeff_tol
) {
  // Same as the R function
  root_info current_root_obj;
  
  Eigen::VectorXd zero_vec = Eigen::VectorXd::Zero(state.size());
  Eigen::VectorXd B = model_spec_obj.lin_root_fun(
    t, 
    zero_vec
  );
  
  Eigen::VectorXd ff3 = model_spec_obj.lin_root_fun(
    t, 
    h * (k1 + k4) - 2.0 * (new_state - state)
  ) - B;
  
  Eigen::VectorXd ff2 = model_spec_obj.lin_root_fun(
    t, 
    3.0 * (new_state - state) - h * (k4 + 2.0 * k1)
  ) - B;
  
  Eigen::VectorXd ff1 = model_spec_obj.lin_root_fun(
    t, 
    k1 * h
  ) - B;
  
  Eigen::VectorXd ff0 = model_spec_obj.lin_root_fun(
    t, 
    state
  );
  
  for (int i = 0; i < ff0.size(); i++) {
    // std::cout << "i: " << i << std::endl;
    
    // if (TRUE) {
    //   NumericVector coeff(4);
    //   coeff(0) = ff0[i];
    //   coeff(1) = ff1[i];
    //   coeff(2) = ff2[i];
    //   coeff(3) = ff3[i];
    //   std::cout << coeff << std::endl;
    // }
    
    // Eigen::VectorXd real_roots = CubicPolyRealRoots( // Cardano's formula, exact, but of course not accurate for ill-conditioned polynomials
    //   ff3[i], ff2[i], ff1[i], ff0[i], precision_real_root_finder
    // );
    // std::cout << std::setprecision(15) << real_roots << std::endl;
    
    // auto roots = boost::math::tools::cubic_roots<double>(ff3[i], ff2[i], ff1[i], ff0[i]); // Using a solver from Boost
    // Eigen::Vector3d real_roots = Eigen::Map<const Eigen::Vector3d>(roots.data());
    
    // Eigen::VectorXd coeffs(4);
    // coeffs << ff0[i], ff1[i], ff2[i], ff3[i];
    
    if (cubic_root_check(ff3[i], ff2[i], ff1[i], ff0[i])) {
      continue; // No roots between 0 and 1 --> can drop the solver step
    }
    
    Eigen::VectorXd coeffs = determine_poly_degree(ff3[i], ff2[i], ff1[i], ff0[i], coeff_tol);
    solver.compute(coeffs);
    Eigen::VectorXcd roots = solver.roots();
    Eigen::Array<bool, Eigen::Dynamic, 1> real_roots_index = 
      roots.imag().array().abs() < precision_real_root_finder;
    Eigen::VectorXd real_roots(real_roots_index.count());
    int k = 0;
    for (int i = 0; i < roots.size(); i++) {
      // std::cout << "t: " << t << std::endl; 
      if (real_roots_index[i]) {
        real_roots[k] = roots[i].real();
        // std::cout << real_roots[k] << std::endl; 
        k += 1; 
      }
    }
    
    // current_root_obj(0) = -1; // which_root
    // current_root_obj(1) = 2.0; // root_time
    // current_root_obj(2) = -1; // root_type
    
    double lower_bound = 0.0; 
    int last_root_root_type = last_root_obj.root_type; 
    int last_root_which_root = last_root_obj.which_root;
    if (last_root_root_type == 1 && last_root_which_root == i) { 
      lower_bound = last_root_offset;
    }
    
    int num_of_allowed_roots = 0;
    for (int j = 0; j < real_roots.size(); j++) {
      if (real_roots[j] >= lower_bound && real_roots[j] <= 1.0) {
        num_of_allowed_roots += 1; 
      }
    }
    
    Eigen::VectorXd allowed_roots(num_of_allowed_roots); 
    int count_allowed_roots = 0; 
    
    if (num_of_allowed_roots > 0) {
      
      for (int j = 0; j < real_roots.size(); j++) {
        if (real_roots[j] >= lower_bound && real_roots[j] <= 1.0) {
          allowed_roots[count_allowed_roots] = real_roots[j];
          count_allowed_roots += 1; 
        }
      }
      
    }
    
    
    // std::cout << allowed_roots << std::endl;
    if (allowed_roots.size() > 0) {
      double min_allowed_roots = allowed_roots.minCoeff();
      if (min_allowed_roots < current_root_obj.root_time) {
        current_root_obj.root_time = min_allowed_roots; 
        current_root_obj.which_root = i;
        current_root_obj.root_type = 1;
      }
    }
  }
  
  return current_root_obj; 
  
}

template<class model_spec>
Rcpp::List deSolverRoot(
    model_spec& model_spec_obj,
    const Eigen::VectorXd& z0,
    const Eigen::VectorXd& times,
    const double& h_max,
    const double& atol,
    const double& rtol,
    const double& last_root_offset_lin_root_finder,
    const double& precision_real_root_lin_root_finder,
    const double& coeff_tol,
    const Eigen::VectorXi& relevant_indices,
    std::string trajectory_info,
    const double& progress_ratio
) {
  
  int D_z0 = z0.size();
  int D_relevant_indices = relevant_indices.size();
  bool equal_dimension_indicator;
  int D;
  if (D_relevant_indices == D_z0) {
    D = D_z0;
    equal_dimension_indicator = TRUE;
  } else {
    D = D_relevant_indices;
    equal_dimension_indicator = FALSE;
  }
  Eigen::VectorXd z = z0;
  double Tmax = times.tail(1)(0);
  double h = h_max;
  int length_times = times.size();
  
  
  //////
  const auto clock_start = std::chrono::steady_clock::now();
  double next_time_to_print = progress_ratio * Tmax;
  //////
  
  Eigen::MatrixXd samples(times.size(), D + 1);
  double t = times(0);
  samples(0, 0) = t;
  if (equal_dimension_indicator) {
    samples.row(0).tail(D) = z.transpose();
  } else {
    for (int k = 1; k < D + 1; k++) {
      // std::cout << "k: " << k << std::endl;
      // std::cout << "relevant_indices_k: " << relevant_indices(k - 1) << std::endl; 
      samples(0, k) = z(relevant_indices(k - 1));
    }
  }
  int samples_count = 1;
  
  Eigen::VectorXd s_times(length_times + 1);
  s_times.segment(0, length_times + 1) = times;
  s_times(length_times) = Tmax + 1000.0;
  
  Eigen::VectorXd k1 = model_spec_obj.ode(
    t,
    z
  );
  
  double PI_alpha = 0.7 / 3.0;
  double PI_beta = 0.4 / 3.0;
  
  root_info last_root_obj;
  
  Eigen::VectorXd z_new; 
  Eigen::VectorXd k4;
  double err; 
  
  int n_evals = 1;
  double old_err = 1.0;
  
  
  while (t < Tmax) {
    bool not_done = TRUE;
    int n_trials = 0;
    h = std::min(h, Tmax - t);
    while (not_done && n_trials < 20) {
      n_trials += 1;
      Eigen::VectorXd ev1 = z + 0.5 * h * k1;
      double t2 = t + 0.5 * h;
      Eigen::VectorXd k2 = model_spec_obj.ode(
        t2,
        ev1
      );
      Eigen::VectorXd ev2 = z + 0.75 * h * k2;
      double t3 = t + 0.75 * h;
      Eigen::VectorXd k3 = model_spec_obj.ode(
        t3,
        ev2
      );
      
      z_new =
        z + h * ((2.0 / 9.0) * k1 + (1.0 / 3.0) * k2 + (4.0 / 9.0) *
        k3);
      double t4 = t + h;
      k4 = model_spec_obj.ode(
        t4,
        z_new
      );
      Eigen::VectorXd z_low =
        z + h * ((7.0 / 24.0) * k1 + 0.25 * k2 + (1.0 / 3.0) * k3 + (1.0 / 8.0) *
        k4);
      err = ((z_new - z_low).array().abs() /
        (atol + rtol * z.array().abs().cwiseMax(z_new.array().abs()))).maxCoeff();
      // std::cout << err << std::endl; 
      n_evals += 3;
      
      if (std::isfinite(err) == 1) {
        
        if (err < 1.0) {
          not_done = FALSE;
          break;
        } else {
          h = h * std::max(0.1, 0.9 * std::pow(err, -1.0 / 3.0));
        }
        
      } else {
        h = 0.1 * h;
      }
    }
    
    if (not_done == TRUE) {
      std::cout << "k1" << std::endl; 
      std::cout << k1 << std::endl;
      std::cout << "z" << std::endl; 
      std::cout << z << std::endl;
      std::cout << "t" << std::endl; 
      std::cout << t << std::endl;
      Rcpp::stop("problem with making progress!"); 
    }
    
    root_info current_root_obj;
    
    root_info rt_out = deSolverLinRootFinder(
      model_spec_obj,
      z,
      z_new,
      k1,
      k4,
      t,
      h,
      last_root_obj,
      last_root_offset_lin_root_finder,
      precision_real_root_lin_root_finder,
      coeff_tol
    );
    
    if (rt_out.root_time < current_root_obj.root_time) {
      current_root_obj = rt_out;
    }
    
    double last = std::min(1.0, current_root_obj.root_time);
    
    // # sets up interpolation polynomial
    Eigen::VectorXd ff3 = h * (k1 + k4) - 2.0 * (z_new - z);
    Eigen::VectorXd ff2 = 3.0 * (z_new - z) - h * (k4 + 2.0 * k1);
    Eigen::VectorXd ff1 = k1 * h;
    
    // collect samples until min(root.time,t+h)
    while (t + h * last >= s_times(samples_count) - 1.0e-14) {
      double r = (s_times(samples_count) - t) / h;
      Eigen::VectorXd samples_row(D + 1);
      samples_row(0) = s_times(samples_count);
      if (equal_dimension_indicator) {
        samples_row.segment(1, D) = z + r * ff1 + std::pow(r, 2) * ff2 + std::pow(r, 3) * ff3;
      } else {
        for (int k = 1; k < D + 1; k++) {
          samples_row(k) = z(relevant_indices(k - 1)) + 
            r * ff1(relevant_indices(k - 1)) + 
            std::pow(r, 2) * ff2(relevant_indices(k - 1)) + 
            std::pow(r, 3) * ff3(relevant_indices(k - 1));
        }
      }
      samples.row(samples_count) = samples_row;
      samples_count += 1;
    }
    
    // do event if required
    
    if (current_root_obj.which_root > -0.1) { // which_root starts from zero in C++
      Eigen::VectorXd z_at_event = z + last * ff1 + std::pow(last, 2) * ff2 + std::pow(last, 3) * ff3;
      double t_new = t + h * last;
      z_new = model_spec_obj.region_lin_event_fun(
        t_new,
        z_at_event
      );
      k4 = model_spec_obj.ode(
        t_new, 
        z_new
      );
      n_evals += 1;
    }
    
    // prepare for next step
    k1 = k4;
    t = t + h * last;
    z = z_new;
    last_root_obj = current_root_obj;
    
    // adapt h using PI controller
    if (err > 1.0e-6) {
      double h_fac = 0.95 * std::pow(err, -PI_alpha) * std::pow(old_err, PI_beta);
      h = h * std::min(2.0, std::max(0.2, h_fac));
    } else {
      h = h_max; // new safeguard to handle linear dynamics
    }
    
    h = std::min(h, h_max);
    // std::cout << h << std::endl; 
    old_err = err;
    
    //// Progress bar ////
    if (t >= next_time_to_print) {
      const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - clock_start).count();
        Rprintf("[%s] t = %.6g / %.6g (%.1f%%), %.0f s elapsed, ~%.0f s left \n",
            trajectory_info.c_str(), t, Tmax, 100.0 * t / Tmax, secs,
                secs * (Tmax - t) / t);
        R_FlushConsole();
        while (next_time_to_print <= t) {
          next_time_to_print += progress_ratio * Tmax;  
        }
    }
    ////
    
  }
  
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("samples") = samples
  );
  output["n_evals"] = n_evals;
  output["n_evals_ode"] = model_spec_obj.n_evals_ode;
  output["n_boundary_crossings"] = model_spec_obj.n_boundary_crossings;
  output["z_final"] = z;
  
  return output;
}

template<class model_list>
Rcpp::List run_sampler( // a general wrapper used for any model defined in C++, see a model example.
    model_list& model_list_obj,
    Rcpp::List settings
) {
  
  Eigen::VectorXd zbar0 = Rcpp::as<Eigen::VectorXd>(settings["zbar0"]);
  int n_parameters = Rcpp::as<int>(settings["n_parameters"]);
  Eigen::VectorXd s_elements = Rcpp::as<Eigen::VectorXd>(settings["s_elements"]);
  Eigen::VectorXd m_elements = Rcpp::as<Eigen::VectorXd>(settings["m_elements"]);
  
  // Set the starting region before the first gradient evaluation
  Eigen::VectorXd qbar0 = zbar0.segment(3, n_parameters);
  model_list_obj.grad_jump_fun(qbar0.cwiseProduct(s_elements) + m_elements);
  
  model_spec_fun<model_list> model_spec_obj( // create a model_spec_obj to define ode, lin_root_fun etc. that can then be used in DeSolverRoot. 
    model_list_obj, 
    settings["lambda"], 
    settings["n_parameters"], 
    settings["push_eps_tol"],
    settings["m_elements"],
    settings["s_elements"]
  );
  
  Rcpp::List output = deSolverRoot(
    model_spec_obj,
    settings["zbar0"],
    settings["times"],
    settings["h_max"],
    settings["atol"],
    settings["rtol"],
    settings["last_root_offset_lin_root_finder"],
    settings["precision_real_root_lin_root_finder"],
    settings["coeff_tol"],
    settings["relevant_indices"],
    settings["trajectory_info"],
    settings["progress_ratio"]
  );
  
  return output;
  
}

}
