#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

//' Build M matrix using C++ (vectorized)
//'
//' Fast C++ implementation of M matrix construction for RC model.
//'
//' @param DeltaD_g NumericVector of DeltaD for group g
//' @param nb_u Integer, number of usable observations
//' @param K Integer, number of lags
//' @return NumericMatrix M_g (nb_u x (K+1))
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix build_M_matrix_cpp(NumericVector DeltaD_g, int nb_u, int K) {
  NumericMatrix M_g(nb_u, K + 1);

  // Original logic: M_g[i, ] <- rev(DeltaD_g[(i + 1):(i + K + 1)])
  // After reversal, column k+1 gets element from position (i + 1 + K - k)
  // C++ uses 0-based indexing, R uses 1-based

  for (int k = 0; k <= K; k++) {
    int start_idx = 1 + K - k;  // R indexing: 2 + K - k, but 0-based in C++
    for (int i = 0; i < nb_u; i++) {
      M_g(i, k) = DeltaD_g[start_idx + i];
    }
  }

  return M_g;
}


//' Compute projection matrix Pi_g efficiently in C++
//'
//' Computes Pi_g = I - M_g %*% ginv(M_g) using Eigen for stability
//'
//' @param M_g NumericMatrix M for group g
//' @return NumericMatrix projection matrix Pi_g
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd compute_projection_cpp(Eigen::Map<Eigen::MatrixXd> M_g) {
  int nb_u = M_g.rows();

  // Use SVD for pseudo-inverse (more stable than ginv)
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(M_g, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Set tolerance for singular values
  double tolerance = 1e-10 * svd.singularValues()(0);

  // Compute pseudo-inverse using SVD
  Eigen::MatrixXd M_pinv = svd.matrixV() *
    (svd.singularValues().array() > tolerance).select(
      svd.singularValues().array().inverse(), 0
    ).matrix().asDiagonal() *
    svd.matrixU().transpose();

  // Compute projection matrix: I - M * M_pinv
  Eigen::MatrixXd Pi_g = Eigen::MatrixXd::Identity(nb_u, nb_u) - M_g * M_pinv;

  return Pi_g;
}


//' Extract group data efficiently in C++
//'
//' Fast extraction of group-specific data
//'
//' @param data NumericVector or NumericMatrix
//' @param indices IntegerVector of indices to extract
//' @return Extracted data
//' @keywords internal
// [[Rcpp::export]]
NumericVector extract_by_indices_cpp(NumericVector data, IntegerVector indices) {
  int n = indices.size();
  NumericVector result(n);

  for (int i = 0; i < n; i++) {
    result[i] = data[indices[i] - 1];  // R uses 1-based indexing
  }

  return result;
}


//' Compute weighted projection and store results (Phase 1 loop body)
//'
//' Optimized C++ implementation of the inner Phase 1 loop
//'
//' @param DeltaY_g NumericVector of DeltaY for group g
//' @param X_g NumericMatrix of X for group g (can be empty)
//' @param DeltaD_g NumericVector of DeltaD for group g
//' @param nb_u Integer, number of usable observations
//' @param K Integer, number of lags
//' @param sqrt_weight Double, square root of weight for this group
//' @return List containing M_g, newY_g, newX_g
//' @keywords internal
// [[Rcpp::export]]
List phase1_group_cpp(NumericVector DeltaY_g, NumericMatrix X_g,
                      NumericVector DeltaD_g, int nb_u, int K, double sqrt_weight) {

  // Build M_g matrix
  NumericMatrix M_g = build_M_matrix_cpp(DeltaD_g, nb_u, K);

  // Convert to Eigen for matrix operations
  Eigen::Map<Eigen::MatrixXd> M_eigen(as<Eigen::Map<Eigen::MatrixXd>>(M_g));
  Eigen::Map<Eigen::VectorXd> Y_eigen(as<Eigen::Map<Eigen::VectorXd>>(DeltaY_g));

  // Compute projection matrix
  Eigen::MatrixXd Pi_g = compute_projection_cpp(M_eigen);

  // Apply projection to Y
  Eigen::VectorXd newY_g = sqrt_weight * (Pi_g * Y_eigen);

  // Apply projection to X if present
  Eigen::MatrixXd newX_g;
  if (X_g.ncol() > 0) {
    Eigen::Map<Eigen::MatrixXd> X_eigen(as<Eigen::Map<Eigen::MatrixXd>>(X_g));
    newX_g = sqrt_weight * (Pi_g * X_eigen);
  }

  return List::create(
    Named("M_g") = M_g,
    Named("newY_g") = newY_g,
    Named("newX_g") = newX_g
  );
}


//' Compute beta coefficients for group (Phase 2 loop body)
//'
//' Optimized C++ implementation of Phase 2 beta computation
//'
//' @param M_g NumericMatrix M for group g
//' @param DeltaY_Xg_gamma NumericVector residuals after removing X*gamma
//' @param K Integer, number of lags
//' @param weight_j Double, weight for this group
//' @return List with ind_avg and beta_hat vectors
//' @keywords internal
// [[Rcpp::export]]
List phase2_group_cpp(NumericMatrix M_g, NumericVector DeltaY_Xg_gamma,
                      int K, double weight_j) {

  int K1 = K + 1;
  NumericVector ind_avg(K1);
  NumericVector beta_hat(K1);

  // Convert to Eigen
  Eigen::Map<Eigen::MatrixXd> M_eigen(as<Eigen::Map<Eigen::MatrixXd>>(M_g));
  Eigen::Map<Eigen::VectorXd> residuals_eigen(as<Eigen::Map<Eigen::VectorXd>>(DeltaY_Xg_gamma));

  // Compute pseudo-inverse of M_g^T
  Eigen::MatrixXd Mt = M_eigen.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mt, Eigen::ComputeThinU | Eigen::ComputeThinV);
  double tolerance = 1e-10 * svd.singularValues()(0);

  Eigen::MatrixXd pinv_Mt = svd.matrixV() *
    (svd.singularValues().array() > tolerance).select(
      svd.singularValues().array().inverse(), 0
    ).matrix().asDiagonal() *
    svd.matrixU().transpose();

  // For each k, compute beta
  for (int k = 0; k <= K; k++) {
    // Create unit vector e_k
    Eigen::VectorXd e_k = Eigen::VectorXd::Zero(K1);
    e_k(k) = 1.0;

    // Compute f_gk = pinv(M^T) * e_k
    Eigen::VectorXd f_gk = pinv_Mt * e_k;

    // Check identifiability: ||M^T * f_gk - e_k|| < 1e-8
    Eigen::VectorXd check = Mt * f_gk - e_k;
    double norm_check = check.norm();

    if (norm_check < 1e-8) {
      ind_avg[k] = weight_j;
      beta_hat[k] = f_gk.dot(residuals_eigen);
    }
  }

  return List::create(
    Named("ind_avg") = ind_avg,
    Named("beta_hat") = beta_hat
  );
}


//' Full RC model estimation in C++ (base model)
//'
//' Complete C++ implementation of the base RC model estimation.
//' This replaces the entire R loop structure with optimized C++ code.
//'
//' @param K Integer, number of lags
//' @param list_g_data IntegerVector, group identifiers for each observation
//' @param list_g IntegerVector, unique group values
//' @param DeltaY NumericVector, outcome differences
//' @param DeltaD NumericVector, treatment differences
//' @param D NumericVector, treatment levels
//' @param X NumericMatrix, covariates (can have 0 columns)
//' @param weights NumericVector, weights for each group
//' @return List with estimation results
//' @keywords internal
// [[Rcpp::export]]
List estim_RC_model_cpp_full(int K, IntegerVector list_g_data, IntegerVector list_g,
                             NumericVector DeltaY, NumericVector DeltaD,
                             NumericVector D, NumericMatrix X, NumericVector weights,
                             bool same_sample = false) {

  int nb_c = list_g.size();
  int n = list_g_data.size();
  int p = X.ncol();

  // Pre-compute group sizes and total observations needed
  IntegerVector group_sizes(nb_c);
  for (int j = 0; j < nb_c; j++) {
    int count = 0;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == list_g[j]) count++;
    }
    group_sizes[j] = count;
  }

  IntegerVector nb_u_vec(nb_c);
  int total_obs = 0;
  for (int j = 0; j < nb_c; j++) {
    nb_u_vec[j] = std::max(0, group_sizes[j] - K - 1);
    total_obs += nb_u_vec[j];
  }

  // Pre-allocate output arrays
  NumericVector newY(total_obs);
  NumericMatrix newX(total_obs, p);
  NumericMatrix M(total_obs, K + 1);

  // PHASE 1: Compute time effects
  int counter = 0;

  for (int j = 0; j < nb_c; j++) {
    int g = list_g[j];

    // Find indices for this group
    std::vector<int> ind_g_vec;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == g) {
        ind_g_vec.push_back(i);
      }
    }

    int nb_obs_g = ind_g_vec.size();
    int nb_u = nb_obs_g - K - 1;

    if (nb_u >= 1) {
      // Extract data for this group
      NumericVector DeltaY_g(nb_u);
      NumericMatrix X_g(nb_u, p);
      NumericVector DeltaD_g(nb_obs_g);

      // Extract DeltaD for full group
      for (int i = 0; i < nb_obs_g; i++) {
        DeltaD_g[i] = DeltaD[ind_g_vec[i]];
      }

      // Extract Y and X for usable observations (after K+1 initial periods)
      for (int i = 0; i < nb_u; i++) {
        int idx = ind_g_vec[K + 1 + i];
        DeltaY_g[i] = DeltaY[idx];
        for (int col = 0; col < p; col++) {
          X_g(i, col) = X(idx, col);
        }
      }

      // Process this group
      List result = phase1_group_cpp(DeltaY_g, X_g, DeltaD_g, nb_u, K, sqrt(weights[j]));

      NumericMatrix M_g = result["M_g"];
      NumericVector newY_g = result["newY_g"];

      // Store results
      for (int i = 0; i < nb_u; i++) {
        newY[counter + i] = newY_g[i];
        for (int k = 0; k <= K; k++) {
          M(counter + i, k) = M_g(i, k);
        }
      }

      if (p > 0) {
        NumericMatrix newX_g = result["newX_g"];
        for (int i = 0; i < nb_u; i++) {
          for (int col = 0; col < p; col++) {
            newX(counter + i, col) = newX_g(i, col);
          }
        }
      }

      counter += nb_u;
    }
  }

  // Compute gamma using crossprod
  NumericVector gamma(p);
  if (p > 0) {
    Eigen::Map<Eigen::MatrixXd> newX_eigen(as<Eigen::Map<Eigen::MatrixXd>>(newX));
    Eigen::Map<Eigen::VectorXd> newY_eigen(as<Eigen::Map<Eigen::VectorXd>>(newY));

    Eigen::MatrixXd XtX = newX_eigen.transpose() * newX_eigen;
    Eigen::VectorXd XtY = newX_eigen.transpose() * newY_eigen;

    Eigen::VectorXd gamma_eigen = XtX.ldlt().solve(XtY);
    gamma = wrap(gamma_eigen);
  }

  // PHASE 2: Compute coefficients
  counter = 0;
  NumericMatrix ind_avg(nb_c, K + 1);
  NumericMatrix beta_bar_hat(nb_c, K + 1);

  for (int j = 0; j < nb_c; j++) {
    int g = list_g[j];

    // Find indices for this group
    std::vector<int> ind_g_vec;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == g) {
        ind_g_vec.push_back(i);
      }
    }

    int nb_obs_g = ind_g_vec.size();
    int nb_u = nb_obs_g - K - 1;

    if (nb_u >= 1) {
      // Extract and compute residuals
      NumericVector DeltaY_Xg_gamma(nb_u);
      for (int i = 0; i < nb_u; i++) {
        int idx = ind_g_vec[K + 1 + i];
        DeltaY_Xg_gamma[i] = DeltaY[idx];
        for (int col = 0; col < p; col++) {
          DeltaY_Xg_gamma[i] -= X(idx, col) * gamma[col];
        }
      }

      // Get M_g from Phase 1
      NumericMatrix M_g(nb_u, K + 1);
      for (int i = 0; i < nb_u; i++) {
        for (int k = 0; k <= K; k++) {
          M_g(i, k) = M(counter + i, k);
        }
      }

      // Compute betas for this group
      List result = phase2_group_cpp(M_g, DeltaY_Xg_gamma, K, weights[j]);
      NumericVector ind_avg_j = result["ind_avg"];
      NumericVector beta_hat_j = result["beta_hat"];

      for (int k = 0; k <= K; k++) {
        ind_avg(j, k) = ind_avg_j[k];
        beta_bar_hat(j, k) = beta_hat_j[k];
      }

      counter += nb_u;
    }
  }

  // Final aggregation
  NumericVector Nobs(K + 1);
  NumericVector beta_bar_hat_final(K + 1);
  int n_coefs = K + 1;

  if (same_sample) {
    std::vector<bool> group_valid(nb_c, true);
    for (int j = 0; j < nb_c; j++)
      for (int k = 0; k < n_coefs; k++)
        if (ind_avg(j, k) == 0) { group_valid[j] = false; break; }

    int n_valid = 0;
    for (int j = 0; j < nb_c; j++) if (group_valid[j]) n_valid++;

    for (int k = 0; k < n_coefs; k++) {
      double sum_ind = 0.0, sum_weighted = 0.0;
      for (int j = 0; j < nb_c; j++) {
        if (group_valid[j]) {
          sum_ind += ind_avg(j, k);
          sum_weighted += ind_avg(j, k) * beta_bar_hat(j, k);
        }
      }
      Nobs[k] = sum_ind;
      if (sum_ind > 0) beta_bar_hat_final[k] = sum_weighted / sum_ind;
    }

    return List::create(
      Named("gamma") = gamma,
      Named("Nobs") = Nobs,
      Named("B_hat") = beta_bar_hat_final,
      Named("n_valid_groups") = n_valid
    );
  } else {
    for (int k = 0; k < n_coefs; k++) {
      double sum_ind = 0.0, sum_weighted = 0.0;
      for (int j = 0; j < nb_c; j++) {
        sum_ind += ind_avg(j, k);
        sum_weighted += ind_avg(j, k) * beta_bar_hat(j, k);
      }
      Nobs[k] = sum_ind;
      if (sum_ind > 0) beta_bar_hat_final[k] = sum_weighted / sum_ind;
    }

    return List::create(
      Named("gamma") = gamma,
      Named("Nobs") = Nobs,
      Named("B_hat") = beta_bar_hat_final
    );
  }
}


//' Build lower triangular M matrix for full dynamics model
//'
//' @param DeltaD_g NumericVector of DeltaD for group g
//' @param T_periods Integer, number of periods
//' @return NumericMatrix M_g (T_periods-1 x T_periods-1)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix build_M_matrix_full_cpp(NumericVector DeltaD_g, int T_periods) {
  int dim = T_periods - 1;
  NumericMatrix M_g(dim, dim);

  // Lower triangular: M_g(i, k) = DeltaD_g[i - k + 1] for k = 0..i
  for (int i = 0; i < dim; i++) {
    for (int k = 0; k <= i; k++) {
      M_g(i, k) = DeltaD_g[i - k + 1];
    }
  }

  return M_g;
}


//' Full dynamics RC model estimation in C++
//'
//' @param T_periods Integer, number of periods
//' @param list_g_data IntegerVector, group identifiers for each observation
//' @param list_g IntegerVector, unique group values
//' @param DeltaY NumericVector, outcome differences
//' @param DeltaD NumericVector, treatment differences
//' @param X NumericMatrix, covariates
//' @param weights NumericVector, weights for each group
//' @param group_sizes IntegerVector, number of observations per group
//' @return List with estimation results
//' @keywords internal
// [[Rcpp::export]]
List estim_RC_model_full_cpp(int T_periods, IntegerVector list_g_data, IntegerVector list_g,
                              NumericVector DeltaY, NumericVector DeltaD,
                              NumericMatrix X, NumericVector weights,
                              IntegerVector group_sizes, bool same_sample = false) {

  int nb_c = list_g.size();
  int n = list_g_data.size();
  int p = X.ncol();
  int dim = T_periods - 1;

  // Count groups with full observations
  int n_full_groups = 0;
  for (int j = 0; j < nb_c; j++) {
    if (group_sizes[j] == T_periods) n_full_groups++;
  }

  int total_obs = n_full_groups * dim;

  // Pre-allocate output arrays
  NumericVector newY(total_obs);
  NumericMatrix newX(total_obs, p);
  NumericMatrix M(total_obs, dim);

  // PHASE 1: Compute time effects
  int counter = 0;

  for (int j = 0; j < nb_c; j++) {
    int g = list_g[j];

    if (group_sizes[j] != T_periods) continue;

    // Find indices for this group
    std::vector<int> ind_g_vec;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == g) {
        ind_g_vec.push_back(i);
      }
    }

    // Extract DeltaD for this group
    NumericVector DeltaD_g(T_periods);
    for (int i = 0; i < T_periods; i++) {
      DeltaD_g[i] = DeltaD[ind_g_vec[i]];
    }

    // Extract Y and X (skip first observation)
    NumericVector DeltaY_g(dim);
    NumericMatrix X_g(dim, p);
    for (int i = 0; i < dim; i++) {
      int idx = ind_g_vec[i + 1];
      DeltaY_g[i] = DeltaY[idx];
      for (int col = 0; col < p; col++) {
        X_g(i, col) = X(idx, col);
      }
    }

    // Build lower triangular M_g
    NumericMatrix M_g = build_M_matrix_full_cpp(DeltaD_g, T_periods);

    // Check if never-switcher (all DeltaD = 0)
    bool is_never_switcher = true;
    for (int i = 0; i < T_periods; i++) {
      if (DeltaD_g[i] != 0) {
        is_never_switcher = false;
        break;
      }
    }

    double sqrt_w = sqrt(weights[j]);

    if (is_never_switcher) {
      // For never-switchers, use data directly (no projection)
      for (int i = 0; i < dim; i++) {
        newY[counter + i] = sqrt_w * DeltaY_g[i];
        for (int col = 0; col < p; col++) {
          newX(counter + i, col) = sqrt_w * X_g(i, col);
        }
        for (int k = 0; k < dim; k++) {
          M(counter + i, k) = M_g(i, k);
        }
      }
    } else {
      // Apply projection
      Eigen::Map<Eigen::MatrixXd> M_eigen(as<Eigen::Map<Eigen::MatrixXd>>(M_g));
      Eigen::MatrixXd Pi_g = compute_projection_cpp(M_eigen);

      Eigen::Map<Eigen::VectorXd> Y_eigen(as<Eigen::Map<Eigen::VectorXd>>(DeltaY_g));
      Eigen::VectorXd newY_g = sqrt_w * (Pi_g * Y_eigen);

      for (int i = 0; i < dim; i++) {
        newY[counter + i] = newY_g(i);
        for (int k = 0; k < dim; k++) {
          M(counter + i, k) = M_g(i, k);
        }
      }

      if (p > 0) {
        Eigen::Map<Eigen::MatrixXd> X_eigen(as<Eigen::Map<Eigen::MatrixXd>>(X_g));
        Eigen::MatrixXd newX_g = sqrt_w * (Pi_g * X_eigen);
        for (int i = 0; i < dim; i++) {
          for (int col = 0; col < p; col++) {
            newX(counter + i, col) = newX_g(i, col);
          }
        }
      }
    }

    counter += dim;
  }

  // Compute gamma
  NumericVector gamma(p);
  if (p > 0 && counter > 0) {
    Eigen::Map<Eigen::MatrixXd> newX_eigen(as<Eigen::Map<Eigen::MatrixXd>>(newX));
    Eigen::Map<Eigen::VectorXd> newY_eigen(as<Eigen::Map<Eigen::VectorXd>>(newY));

    Eigen::MatrixXd XtX = newX_eigen.transpose() * newX_eigen;
    Eigen::VectorXd XtY = newX_eigen.transpose() * newY_eigen;

    Eigen::VectorXd gamma_eigen = XtX.ldlt().solve(XtY);
    gamma = wrap(gamma_eigen);
  }

  // PHASE 2: Compute coefficients
  counter = 0;
  NumericMatrix ind_avg(nb_c, dim);
  NumericMatrix beta_bar_hat(nb_c, dim);

  for (int j = 0; j < nb_c; j++) {
    int g = list_g[j];

    if (group_sizes[j] != T_periods) continue;

    // Find indices for this group
    std::vector<int> ind_g_vec;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == g) {
        ind_g_vec.push_back(i);
      }
    }

    // Extract DeltaD to check never-switcher
    bool is_never_switcher = true;
    for (int i = 0; i < T_periods; i++) {
      if (DeltaD[ind_g_vec[i]] != 0) {
        is_never_switcher = false;
        break;
      }
    }

    // Compute residuals
    NumericVector DeltaY_Xg_gamma(dim);
    for (int i = 0; i < dim; i++) {
      int idx = ind_g_vec[i + 1];
      DeltaY_Xg_gamma[i] = DeltaY[idx];
      for (int col = 0; col < p; col++) {
        DeltaY_Xg_gamma[i] -= X(idx, col) * gamma[col];
      }
    }

    if (is_never_switcher) {
      // For never-switchers, all coefficients identified with beta = 0
      for (int k = 0; k < dim; k++) {
        ind_avg(j, k) = weights[j];
        beta_bar_hat(j, k) = 0.0;
      }
    } else {
      // Get M_g from Phase 1
      NumericMatrix M_g(dim, dim);
      for (int i = 0; i < dim; i++) {
        for (int k = 0; k < dim; k++) {
          M_g(i, k) = M(counter + i, k);
        }
      }

      // Compute betas
      Eigen::Map<Eigen::MatrixXd> M_eigen(as<Eigen::Map<Eigen::MatrixXd>>(M_g));
      Eigen::Map<Eigen::VectorXd> residuals_eigen(as<Eigen::Map<Eigen::VectorXd>>(DeltaY_Xg_gamma));

      Eigen::MatrixXd Mt = M_eigen.transpose();
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mt, Eigen::ComputeThinU | Eigen::ComputeThinV);
      double tolerance = 1e-10 * svd.singularValues()(0);

      Eigen::MatrixXd pinv_Mt = svd.matrixV() *
        (svd.singularValues().array() > tolerance).select(
          svd.singularValues().array().inverse(), 0
        ).matrix().asDiagonal() *
        svd.matrixU().transpose();

      for (int k = 0; k < dim; k++) {
        Eigen::VectorXd e_k = Eigen::VectorXd::Zero(dim);
        e_k(k) = 1.0;

        Eigen::VectorXd f_gk = pinv_Mt * e_k;
        Eigen::VectorXd check = Mt * f_gk - e_k;

        if (check.norm() < 1e-8) {
          ind_avg(j, k) = weights[j];
          beta_bar_hat(j, k) = f_gk.dot(residuals_eigen);
        }
      }
    }

    counter += dim;
  }

  // Final aggregation
  NumericVector Nobs(dim);
  NumericVector beta_bar_hat_final(dim);

  if (same_sample) {
    std::vector<bool> group_valid(nb_c, true);
    for (int j = 0; j < nb_c; j++)
      for (int k = 0; k < dim; k++)
        if (ind_avg(j, k) == 0) { group_valid[j] = false; break; }

    int n_valid = 0;
    for (int j = 0; j < nb_c; j++) if (group_valid[j]) n_valid++;

    for (int k = 0; k < dim; k++) {
      double sum_ind = 0.0, sum_weighted = 0.0;
      for (int j = 0; j < nb_c; j++) {
        if (group_valid[j]) {
          sum_ind += ind_avg(j, k);
          sum_weighted += ind_avg(j, k) * beta_bar_hat(j, k);
        }
      }
      Nobs[k] = sum_ind;
      if (sum_ind > 0) beta_bar_hat_final[k] = sum_weighted / sum_ind;
    }

    return List::create(
      Named("gamma") = gamma,
      Named("Nobs") = Nobs,
      Named("B_hat") = beta_bar_hat_final,
      Named("n_valid_groups") = n_valid
    );
  } else {
    for (int k = 0; k < dim; k++) {
      double sum_ind = 0.0, sum_weighted = 0.0;
      for (int j = 0; j < nb_c; j++) {
        sum_ind += ind_avg(j, k);
        sum_weighted += ind_avg(j, k) * beta_bar_hat(j, k);
      }
      Nobs[k] = sum_ind;
      if (sum_ind > 0) beta_bar_hat_final[k] = sum_weighted / sum_ind;
    }

    return List::create(
      Named("gamma") = gamma,
      Named("Nobs") = Nobs,
      Named("B_hat") = beta_bar_hat_final
    );
  }
}


//' Build M matrix with interactions for interactions model
//'
//' @param DeltaD_g NumericVector of DeltaD for group g
//' @param D_g NumericVector of D for group g
//' @param nb_u Integer, number of usable observations
//' @param K Integer, number of lags
//' @return NumericMatrix M_g with main effects and interactions
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix build_M_matrix_interactions_cpp(NumericVector DeltaD_g, NumericVector D_g,
                                               int nb_u, int K) {
  int n_params = (K + 1) * (K + 2) / 2;
  NumericMatrix M_g(nb_u, n_params);

  int col_idx = 0;

  // Main effects: columns 0 to K
  for (int k = 0; k <= K; k++) {
    for (int i = 0; i < nb_u; i++) {
      M_g(i, col_idx) = DeltaD_g[i + K - k + 1];
    }
    col_idx++;
  }

  // Interaction effects: for k < kp
  for (int k = 0; k < K; k++) {
    for (int kp = k + 1; kp <= K; kp++) {
      for (int i = 0; i < nb_u; i++) {
        // Delta(D_{t-k} * D_{t-k'})
        double curr_prod = D_g[i + K + 1 - k] * D_g[i + K + 1 - kp];
        double prev_prod = D_g[i + K - k] * D_g[i + K - kp];
        M_g(i, col_idx) = curr_prod - prev_prod;
      }
      col_idx++;
    }
  }

  return M_g;
}


//' Interactions RC model estimation in C++
//'
//' @param K Integer, number of lags
//' @param list_g_data IntegerVector, group identifiers for each observation
//' @param list_g IntegerVector, unique group values
//' @param DeltaY NumericVector, outcome differences
//' @param DeltaD NumericVector, treatment differences
//' @param D NumericVector, treatment levels
//' @param X NumericMatrix, covariates
//' @param weights NumericVector, weights for each group
//' @return List with estimation results
//' @keywords internal
// [[Rcpp::export]]
List estim_RC_model_interactions_cpp(int K, IntegerVector list_g_data, IntegerVector list_g,
                                      NumericVector DeltaY, NumericVector DeltaD,
                                      NumericVector D, NumericMatrix X, NumericVector weights,
                                      bool same_sample = false) {

  int nb_c = list_g.size();
  int n = list_g_data.size();
  int p = X.ncol();
  int n_params = (K + 1) * (K + 2) / 2;

  // Pre-compute group sizes
  IntegerVector group_sizes(nb_c);
  for (int j = 0; j < nb_c; j++) {
    int count = 0;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == list_g[j]) count++;
    }
    group_sizes[j] = count;
  }

  // Compute total observations
  int total_obs = 0;
  for (int j = 0; j < nb_c; j++) {
    int nb_u = std::max(0, group_sizes[j] - K - 1);
    total_obs += nb_u;
  }

  // Pre-allocate output arrays
  NumericVector newY(total_obs);
  NumericMatrix newX(total_obs, p);
  NumericMatrix M(total_obs, n_params);

  // PHASE 1: Compute time effects
  int counter = 0;

  for (int j = 0; j < nb_c; j++) {
    int g = list_g[j];
    int nb_obs_g = group_sizes[j];
    int nb_u = nb_obs_g - K - 1;

    if (nb_u < 1) continue;

    // Find indices for this group
    std::vector<int> ind_g_vec;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == g) {
        ind_g_vec.push_back(i);
      }
    }

    // Extract data for this group
    NumericVector DeltaD_g(nb_obs_g);
    NumericVector D_g(nb_obs_g);
    for (int i = 0; i < nb_obs_g; i++) {
      DeltaD_g[i] = DeltaD[ind_g_vec[i]];
      D_g[i] = D[ind_g_vec[i]];
    }

    NumericVector DeltaY_g(nb_u);
    NumericMatrix X_g(nb_u, p);
    for (int i = 0; i < nb_u; i++) {
      int idx = ind_g_vec[K + 1 + i];
      DeltaY_g[i] = DeltaY[idx];
      for (int col = 0; col < p; col++) {
        X_g(i, col) = X(idx, col);
      }
    }

    // Build M_g with interactions
    NumericMatrix M_g = build_M_matrix_interactions_cpp(DeltaD_g, D_g, nb_u, K);

    // Apply projection
    Eigen::Map<Eigen::MatrixXd> M_eigen(as<Eigen::Map<Eigen::MatrixXd>>(M_g));
    Eigen::MatrixXd Pi_g = compute_projection_cpp(M_eigen);

    double sqrt_w = sqrt(weights[j]);

    Eigen::Map<Eigen::VectorXd> Y_eigen(as<Eigen::Map<Eigen::VectorXd>>(DeltaY_g));
    Eigen::VectorXd newY_g = sqrt_w * (Pi_g * Y_eigen);

    for (int i = 0; i < nb_u; i++) {
      newY[counter + i] = newY_g(i);
      for (int k = 0; k < n_params; k++) {
        M(counter + i, k) = M_g(i, k);
      }
    }

    if (p > 0) {
      Eigen::Map<Eigen::MatrixXd> X_eigen(as<Eigen::Map<Eigen::MatrixXd>>(X_g));
      Eigen::MatrixXd newX_g = sqrt_w * (Pi_g * X_eigen);
      for (int i = 0; i < nb_u; i++) {
        for (int col = 0; col < p; col++) {
          newX(counter + i, col) = newX_g(i, col);
        }
      }
    }

    counter += nb_u;
  }

  // Compute gamma
  NumericVector gamma(p);
  if (p > 0 && counter > 0) {
    Eigen::Map<Eigen::MatrixXd> newX_eigen(as<Eigen::Map<Eigen::MatrixXd>>(newX));
    Eigen::Map<Eigen::VectorXd> newY_eigen(as<Eigen::Map<Eigen::VectorXd>>(newY));

    Eigen::MatrixXd XtX = newX_eigen.transpose() * newX_eigen;
    Eigen::VectorXd XtY = newX_eigen.transpose() * newY_eigen;

    Eigen::VectorXd gamma_eigen = XtX.ldlt().solve(XtY);
    gamma = wrap(gamma_eigen);
  }

  // PHASE 2: Compute coefficients
  counter = 0;
  NumericMatrix ind_avg(nb_c, n_params);
  NumericMatrix beta_bar_hat(nb_c, n_params);

  for (int j = 0; j < nb_c; j++) {
    int g = list_g[j];
    int nb_obs_g = group_sizes[j];
    int nb_u = nb_obs_g - K - 1;

    if (nb_u < 1) continue;

    // Find indices for this group
    std::vector<int> ind_g_vec;
    for (int i = 0; i < n; i++) {
      if (list_g_data[i] == g) {
        ind_g_vec.push_back(i);
      }
    }

    // Compute residuals
    NumericVector DeltaY_Xg_gamma(nb_u);
    for (int i = 0; i < nb_u; i++) {
      int idx = ind_g_vec[K + 1 + i];
      DeltaY_Xg_gamma[i] = DeltaY[idx];
      for (int col = 0; col < p; col++) {
        DeltaY_Xg_gamma[i] -= X(idx, col) * gamma[col];
      }
    }

    // Get M_g from Phase 1
    NumericMatrix M_g(nb_u, n_params);
    for (int i = 0; i < nb_u; i++) {
      for (int k = 0; k < n_params; k++) {
        M_g(i, k) = M(counter + i, k);
      }
    }

    // Compute betas
    Eigen::Map<Eigen::MatrixXd> M_eigen(as<Eigen::Map<Eigen::MatrixXd>>(M_g));
    Eigen::Map<Eigen::VectorXd> residuals_eigen(as<Eigen::Map<Eigen::VectorXd>>(DeltaY_Xg_gamma));

    Eigen::MatrixXd Mt = M_eigen.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Mt, Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = 1e-10 * svd.singularValues()(0);

    Eigen::MatrixXd pinv_Mt = svd.matrixV() *
      (svd.singularValues().array() > tolerance).select(
        svd.singularValues().array().inverse(), 0
      ).matrix().asDiagonal() *
      svd.matrixU().transpose();

    for (int param_idx = 0; param_idx < n_params; param_idx++) {
      Eigen::VectorXd e_k = Eigen::VectorXd::Zero(n_params);
      e_k(param_idx) = 1.0;

      Eigen::VectorXd f_gk = pinv_Mt * e_k;
      Eigen::VectorXd check = Mt * f_gk - e_k;

      if (check.norm() < 1e-8) {
        ind_avg(j, param_idx) = weights[j];
        beta_bar_hat(j, param_idx) = f_gk.dot(residuals_eigen);
      }
    }

    counter += nb_u;
  }

  // Final aggregation
  NumericVector Nobs(n_params);
  NumericVector beta_bar_hat_final(n_params);

  if (same_sample) {
    std::vector<bool> group_valid(nb_c, true);
    for (int j = 0; j < nb_c; j++)
      for (int k = 0; k < n_params; k++)
        if (ind_avg(j, k) == 0) { group_valid[j] = false; break; }

    int n_valid = 0;
    for (int j = 0; j < nb_c; j++) if (group_valid[j]) n_valid++;

    for (int k = 0; k < n_params; k++) {
      double sum_ind = 0.0, sum_weighted = 0.0;
      for (int j = 0; j < nb_c; j++) {
        if (group_valid[j]) {
          sum_ind += ind_avg(j, k);
          sum_weighted += ind_avg(j, k) * beta_bar_hat(j, k);
        }
      }
      Nobs[k] = sum_ind;
      if (sum_ind > 0) beta_bar_hat_final[k] = sum_weighted / sum_ind;
    }

    return List::create(
      Named("gamma") = gamma,
      Named("Nobs") = Nobs,
      Named("B_hat") = beta_bar_hat_final,
      Named("n_valid_groups") = n_valid
    );
  } else {
    for (int k = 0; k < n_params; k++) {
      double sum_ind = 0.0, sum_weighted = 0.0;
      for (int j = 0; j < nb_c; j++) {
        sum_ind += ind_avg(j, k);
        sum_weighted += ind_avg(j, k) * beta_bar_hat(j, k);
      }
      Nobs[k] = sum_ind;
      if (sum_ind > 0) beta_bar_hat_final[k] = sum_weighted / sum_ind;
    }

    return List::create(
      Named("gamma") = gamma,
      Named("Nobs") = Nobs,
      Named("B_hat") = beta_bar_hat_final
    );
  }
}
