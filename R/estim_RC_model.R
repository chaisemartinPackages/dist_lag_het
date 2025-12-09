#' Unified RC Model Estimation with Multiple Specifications
#'
#' Main wrapper function that allows choosing between different model specifications
#' for estimating treatment effects in complex designs under parallel-trends assumptions.
#' Optionally computes bootstrap standard errors and confidence intervals.
#'
#' @param K Integer. Number of lags (required for "base" and "interactions" models,
#'   ignored for "full_dynamics"). Default is NULL.
#' @param data Data frame containing the panel data.
#' @param group_col Character. Name of the column identifying groups (e.g., counties, individuals).
#' @param deltaY_col Character. Name of the column for ΔY (first difference of outcome).
#' @param deltaD_col Character. Name of the column for ΔD (first difference of treatment).
#' @param D_col Character. Name of the column for D (treatment level).
#' @param X_cols Character vector. Names of columns for covariates.
#' @param weights Numeric vector. Optional weights for estimation (default: uniform weights).
#' @param model Character. One of "base", "full_dynamics", or "interactions". Default is "base".
#' @param bootstrap Logical. If TRUE, compute bootstrap standard errors (default: FALSE).
#' @param B Integer. Number of bootstrap iterations (default: 1000, ignored if bootstrap = FALSE).
#' @param conf_level Numeric. Confidence level for bootstrap intervals (default: 0.95, ignored if bootstrap = FALSE).
#' @param verbose Logical. Show bootstrap progress bar (default: TRUE, ignored if bootstrap = FALSE).
#'
#' @return A list containing:
#'   \item{gamma}{Coefficients for covariates}
#'   \item{B_hat}{Estimated treatment effect coefficients (betas)}
#'   \item{Nobs}{Number of observations used for each coefficient}
#'   \item{model}{Model type used}
#'   If bootstrap = TRUE, also includes:
#'   \item{se}{Bootstrap standard errors}
#'   \item{t_stat}{t-statistics (B_hat / se)}
#'   \item{ci_lower}{Lower bounds of confidence intervals}
#'   \item{ci_upper}{Upper bounds of confidence intervals}
#'   \item{boot_iterations}{Number of successful bootstrap iterations}
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' data <- generate_did_data(n_groups = 50, n_periods = 5)
#'
#' # Estimate base model without bootstrap
#' result <- estim_RC_model_unified(
#'   K = 2,
#'   data = data,
#'   group_col = "group",
#'   deltaY_col = "delta_Y",
#'   deltaD_col = "delta_D",
#'   D_col = "D",
#'   X_cols = c("X1", "X2"),
#'   model = "base"
#' )
#'
#' # Estimate with bootstrap standard errors
#' result_boot <- estim_RC_model_unified(
#'   K = 2,
#'   data = data,
#'   group_col = "group",
#'   deltaY_col = "delta_Y",
#'   deltaD_col = "delta_D",
#'   D_col = "D",
#'   X_cols = c("X1", "X2"),
#'   model = "base",
#'   bootstrap = TRUE,
#'   B = 1000
#' )
#'
#' print(result_boot$B_hat)
#' print(result_boot$se)
#' }
#'
#' @export
#' @importFrom MASS ginv
estim_RC_model_unified <- function(K = NULL, data, group_col, deltaY_col, deltaD_col,
                                   D_col, X_cols, weights = NULL,
                                   model = c("base", "full_dynamics", "interactions"),
                                   bootstrap = FALSE, B = 1000, conf_level = 0.95,
                                   verbose = TRUE) {

  model <- match.arg(model)

  # Estimate the main model
  if (model == "base") {
    if (is.null(K)) stop("K must be specified for base model")
    result <- estim_RC_model_base(K, data, group_col, deltaY_col, deltaD_col,
                               D_col, X_cols, weights)

  } else if (model == "full_dynamics") {
    if (!is.null(K)) warning("K is ignored for full_dynamics model")
    result <- estim_RC_model_full(K = NULL, data, group_col, deltaY_col, deltaD_col,
                               D_col, X_cols, weights)

  } else if (model == "interactions") {
    if (is.null(K)) stop("K must be specified for interactions model")
    result <- estim_RC_model_interactions(K, data, group_col, deltaY_col, deltaD_col,
                                       D_col, X_cols, weights)
  }

  # If bootstrap is requested, compute standard errors
  if (bootstrap) {
    boot_stats <- .bootstrap_RC_model_internal(
      result = result,
      data = data,
      group_col = group_col,
      deltaY_col = deltaY_col,
      deltaD_col = deltaD_col,
      D_col = D_col,
      X_cols = X_cols,
      K = K,
      B = B,
      conf_level = conf_level,
      verbose = verbose
    )

    # Add bootstrap results to output
    result$se <- boot_stats$se
    result$t_stat <- boot_stats$t_stat
    result$ci_lower <- boot_stats$ci_lower
    result$ci_upper <- boot_stats$ci_upper
    result$boot_iterations <- boot_stats$boot_iterations
  }

  return(result)
}


#' Base RC Model Estimation
#'
#' Estimates the base RC model with K lags without interaction terms.
#'
#' @inheritParams estim_RC_model_unified
#'
#' @return A list containing estimation results.
#'
#' @export
#' @importFrom MASS ginv
estim_RC_model_base <- function(K, data, group_col, deltaY_col, deltaD_col,
                               D_col, X_cols, weights = NULL) {

  # Convert to dataframe if needed and ensure numeric columns
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Extract columns by name
  list_g_data <- data[[group_col]]
  list_g <- unique(list_g_data)
  nb_c <- length(list_g)

  # If no weights provided, use uniform weights
  if (is.null(weights)) {
    weights <- rep(1, nb_c)
  }

  # Create combined matrix of DeltaY and X variables
  DeltaY_X <- as.matrix(cbind(data[[deltaY_col]], data[, X_cols, drop = FALSE]))
  DeltaD <- data[[deltaD_col]]
  D <- data[[D_col]]

  # Ensure numeric
  DeltaY_X <- apply(DeltaY_X, 2, as.numeric)
  DeltaD <- as.numeric(DeltaD)
  D <- as.numeric(D)

  ## PHASE 1: Estimation of time effects
  n <- nrow(data)
  p <- length(X_cols)  # Number of X variables
  newY <- numeric(n)
  newX <- matrix(0, n, p)
  M <- matrix(0, n, K + 1)
  j <- 1
  counter <- 0

  for (g in list_g) {
    ind_g <- (list_g_data == g)
    nb_obs_g <- sum(ind_g)
    nb_u <- nb_obs_g - K - 1

    if (nb_u >= 1) {
      DeltaY_X_g <- DeltaY_X[ind_g, , drop = FALSE]
      DeltaY_g <- DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 1]
      X_g <- DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 2:ncol(DeltaY_X_g), drop = FALSE]

      DeltaD_g <- DeltaD[ind_g]
      M_g <- matrix(0, nb_u, K + 1)
      for (i in 1:nb_u) {
        M_g[i, ] <- rev(DeltaD_g[(i + 1):(i + K + 1)])
      }

      Pi_g <- diag(nb_u) - M_g %*% MASS::ginv(M_g)
      M[(counter + 1):(counter + nb_u), ] <- M_g
      newY[(counter + 1):(counter + nb_u)] <- sqrt(weights[j]) * (Pi_g %*% DeltaY_g)
      newX[(counter + 1):(counter + nb_u), ] <- sqrt(weights[j]) * (Pi_g %*% X_g)
      counter <- counter + nb_u
    }
    j <- j + 1
  }

  newY <- newY[1:counter]
  newX <- as.matrix(newX[1:counter, , drop = FALSE])
  M <- as.matrix(M[1:counter, , drop = FALSE])
  gamma <- solve(t(newX) %*% newX, t(newX) %*% newY)

  ## PHASE 2: Estimation of coefficients
  counter <- 0
  j <- 1
  ind_avg <- matrix(0, nb_c, K + 1)
  beta_bar_hat <- matrix(0, nb_c, K + 1)

  for (g in list_g) {
    ind_g <- (list_g_data == g)
    nb_obs_g <- sum(ind_g)
    nb_u <- nb_obs_g - K - 1

    if (nb_u >= 1) {
      DeltaY_X_g <- DeltaY_X[ind_g, , drop = FALSE]
      DeltaY_Xg_gamma <- DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 1] -
        DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 2:ncol(DeltaY_X_g), drop = FALSE] %*% gamma

      M_g <- M[(counter + 1):(counter + nb_u), , drop = FALSE]
      pinv_Mg <- MASS::ginv(t(M_g))

      for (k in 0:K) {
        ek <- numeric(K + 1)
        ek[k + 1] <- 1
        f_gk <- pinv_Mg %*% ek

        if (sqrt(sum((t(M_g) %*% f_gk - ek)^2)) < 1e-8) {
          ind_avg[j, k + 1] <- weights[j]
          beta_bar_hat[j, k + 1] <- as.numeric(t(f_gk) %*% DeltaY_Xg_gamma)
        }
      }
      counter <- counter + max(0, nb_obs_g - K - 1)
    }
    j <- j + 1
  }

  Nobs <- colSums(ind_avg)
  beta_bar_hat_final <- numeric(K + 1)
  for (k in 1:(K + 1)) {
    if (Nobs[k] > 0) {
      beta_bar_hat_final[k] <- sum(ind_avg[, k] * beta_bar_hat[, k]) / Nobs[k]
    }
  }

  names(beta_bar_hat_final) <- paste0("beta_", 0:K)
  return(list(gamma = gamma, Nobs = Nobs, B_hat = beta_bar_hat_final, model = "base"))
}


#' Full Dynamics RC Model Estimation
#'
#' Estimates the RC model with full dynamics (all available periods).
#'
#' @inheritParams estim_RC_model_unified
#'
#' @return A list containing estimation results.
#'
#' @export
#' @importFrom MASS ginv
estim_RC_model_full <- function(K = NULL, data, group_col, deltaY_col, deltaD_col,
                                D_col, X_cols, weights = NULL) {

  # Convert to dataframe if needed and ensure numeric columns
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Extract columns by name
  list_g_data <- data[[group_col]]
  list_g <- unique(list_g_data)
  nb_c <- length(list_g)

  # If no weights provided, use uniform weights
  if (is.null(weights)) {
    weights <- rep(1, nb_c)
  }

  # Create combined matrix of DeltaY and X variables
  DeltaY_X <- as.matrix(cbind(data[[deltaY_col]], data[, X_cols, drop = FALSE]))
  DeltaD <- data[[deltaD_col]]
  D <- data[[D_col]]

  # Ensure numeric
  DeltaY_X <- apply(DeltaY_X, 2, as.numeric)
  DeltaD <- as.numeric(DeltaD)
  D <- as.numeric(D)

  # Determine T (number of periods)
  T_periods <- max(table(list_g_data))

  ## PHASE 1: Estimation of time effects
  n <- nrow(data)
  p <- ncol(DeltaY_X) - 1
  nb_c <- length(list_g)
  newY <- numeric(n)
  newX <- matrix(0, n, p)
  M <- matrix(0, n, T_periods - 1)
  j <- 1
  counter <- 0

  for (g in list_g) {
    ind_g <- (list_g_data == g)
    nb_obs_g <- sum(ind_g)

    if (nb_obs_g == T_periods) {
      # Extract data for this group
      DeltaY_X_g <- DeltaY_X[ind_g, , drop = FALSE]
      DeltaD_g <- DeltaD[ind_g]

      # Skip first observation (we need it for constructing M)
      DeltaY_g <- DeltaY_X_g[2:T_periods, 1]
      X_g <- DeltaY_X_g[2:T_periods, 2:ncol(DeltaY_X_g), drop = FALSE]

      # Construct M_g: lower triangular matrix
      M_g <- matrix(0, T_periods - 1, T_periods - 1)
      for (i in 1:(T_periods - 1)) {
        for (k in 0:(i - 1)) {
          M_g[i, k + 1] <- DeltaD_g[i - k + 1]
        }
      }

      # Check if group is a never-switcher
      if (all(DeltaD_g == 0)) {
        # For never-switchers, use them directly
        newY[(counter + 1):(counter + T_periods - 1)] <- sqrt(weights[j]) * DeltaY_g
        newX[(counter + 1):(counter + T_periods - 1), ] <- sqrt(weights[j]) * X_g
        M[(counter + 1):(counter + T_periods - 1), ] <- M_g
      } else {
        # For switchers, apply projection
        Pi_g <- diag(T_periods - 1) - M_g %*% MASS::ginv(M_g)
        newY[(counter + 1):(counter + T_periods - 1)] <- sqrt(weights[j]) * (Pi_g %*% DeltaY_g)
        newX[(counter + 1):(counter + T_periods - 1), ] <- sqrt(weights[j]) * (Pi_g %*% X_g)
        M[(counter + 1):(counter + T_periods - 1), ] <- M_g
      }
      counter <- counter + T_periods - 1
    }
    j <- j + 1
  }

  newY <- newY[1:counter]
  newX <- as.matrix(newX[1:counter, , drop = FALSE])
  M <- as.matrix(M[1:counter, , drop = FALSE])
  gamma <- solve(t(newX) %*% newX, t(newX) %*% newY)

  ## PHASE 2: Estimation of coefficients
  counter <- 0
  j <- 1
  ind_avg <- matrix(0, nb_c, T_periods - 1)
  beta_bar_hat <- matrix(0, nb_c, T_periods - 1)

  for (g in list_g) {
    ind_g <- (list_g_data == g)
    nb_obs_g <- sum(ind_g)

    if (nb_obs_g == T_periods) {
      DeltaY_X_g <- DeltaY_X[ind_g, , drop = FALSE]
      DeltaY_Xg_gamma <- DeltaY_X_g[2:T_periods, 1] -
        DeltaY_X_g[2:T_periods, 2:ncol(DeltaY_X_g), drop = FALSE] %*% gamma

      M_g <- M[(counter + 1):(counter + T_periods - 1), , drop = FALSE]

      # Check if it's a never-switcher
      if (all(DeltaD[ind_g] == 0)) {
        # For never-switchers, all coefficients are identified
        for (k in 0:(T_periods - 2)) {
          ind_avg[j, k + 1] <- weights[j]
          beta_bar_hat[j, k + 1] <- 0  # Effect is zero for never-switchers
        }
      } else {
        pinv_Mg <- MASS::ginv(t(M_g))

        for (k in 0:(T_periods - 2)) {
          ek <- numeric(T_periods - 1)
          ek[k + 1] <- 1
          f_gk <- pinv_Mg %*% ek

          if (sqrt(sum((t(M_g) %*% f_gk - ek)^2)) < 1e-8) {
            ind_avg[j, k + 1] <- weights[j]
            beta_bar_hat[j, k + 1] <- as.numeric(t(f_gk) %*% DeltaY_Xg_gamma)
          }
        }
      }
      counter <- counter + T_periods - 1
    }
    j <- j + 1
  }

  Nobs <- colSums(ind_avg)
  beta_bar_hat_final <- numeric(T_periods - 1)
  for (k in 1:(T_periods - 1)) {
    if (Nobs[k] > 0) {
      beta_bar_hat_final[k] <- sum(ind_avg[, k] * beta_bar_hat[, k]) / Nobs[k]
    }
  }

  names(beta_bar_hat_final) <- paste0("beta_", 0:(T_periods - 2))
  return(list(gamma = gamma, Nobs = Nobs, B_hat = beta_bar_hat_final, model = "full_dynamics"))
}


#' RC Model with Interaction Terms
#'
#' Estimates the RC model with interaction terms between lags.
#'
#' @inheritParams estim_RC_model_unified
#'
#' @return A list containing estimation results.
#'
#' @export
#' @importFrom MASS ginv
estim_RC_model_interactions <- function(K, data, group_col, deltaY_col, deltaD_col,
                                        D_col, X_cols, weights = NULL) {

  # Convert to dataframe if needed and ensure numeric columns
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Extract columns by name
  list_g_data <- data[[group_col]]
  list_g <- unique(list_g_data)
  nb_c <- length(list_g)

  # If no weights provided, use uniform weights
  if (is.null(weights)) {
    weights <- rep(1, nb_c)
  }

  # Create combined matrix of DeltaY and X variables
  DeltaY_X <- as.matrix(cbind(data[[deltaY_col]], data[, X_cols, drop = FALSE]))
  DeltaD <- data[[deltaD_col]]
  D <- data[[D_col]]

  # Ensure numeric
  DeltaY_X <- apply(DeltaY_X, 2, as.numeric)
  DeltaD <- as.numeric(DeltaD)
  D <- as.numeric(D)

  # Calculate number of parameters: main effects + interactions
  n_params <- (K + 1) * (K + 2) / 2

  ## PHASE 1: Estimation of time effects
  n <- nrow(data)
  p <- ncol(DeltaY_X) - 1
  nb_c <- length(list_g)
  newY <- numeric(n)
  newX <- matrix(0, n, p)
  M <- matrix(0, n, n_params)
  j <- 1
  counter <- 0

  for (g in list_g) {
    ind_g <- (list_g_data == g)
    nb_obs_g <- sum(ind_g)
    nb_u <- nb_obs_g - K - 1

    if (nb_u >= 1) {
      DeltaY_X_g <- DeltaY_X[ind_g, , drop = FALSE]
      DeltaY_g <- DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 1]
      X_g <- DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 2:ncol(DeltaY_X_g), drop = FALSE]

      DeltaD_g <- DeltaD[ind_g]
      D_g <- D[ind_g]

      # Construct M_g with interaction terms
      M_g <- matrix(0, nb_u, n_params)
      col_idx <- 1

      # Main effects
      for (k in 0:K) {
        for (i in 1:nb_u) {
          M_g[i, col_idx] <- DeltaD_g[i + K - k + 1]
        }
        col_idx <- col_idx + 1
      }

      # Interaction effects (for k < k')
      for (k in 0:(K - 1)) {
        for (kp in (k + 1):K) {
          for (i in 1:nb_u) {
            # Delta(D_{t-k} * D_{t-k'})
            curr_prod <- D_g[i + K + 1 - k] * D_g[i + K + 1 - kp]
            prev_prod <- D_g[i + K - k] * D_g[i + K - kp]
            M_g[i, col_idx] <- curr_prod - prev_prod
          }
          col_idx <- col_idx + 1
        }
      }

      Pi_g <- diag(nb_u) - M_g %*% MASS::ginv(M_g)
      M[(counter + 1):(counter + nb_u), ] <- M_g
      newY[(counter + 1):(counter + nb_u)] <- sqrt(weights[j]) * (Pi_g %*% DeltaY_g)
      newX[(counter + 1):(counter + nb_u), ] <- sqrt(weights[j]) * (Pi_g %*% X_g)
      counter <- counter + nb_u
    }
    j <- j + 1
  }

  newY <- newY[1:counter]
  newX <- as.matrix(newX[1:counter, , drop = FALSE])
  M <- as.matrix(M[1:counter, , drop = FALSE])
  gamma <- solve(t(newX) %*% newX, t(newX) %*% newY)

  ## PHASE 2: Estimation of coefficients
  counter <- 0
  j <- 1
  ind_avg <- matrix(0, nb_c, n_params)
  beta_bar_hat <- matrix(0, nb_c, n_params)

  for (g in list_g) {
    ind_g <- (list_g_data == g)
    nb_obs_g <- sum(ind_g)
    nb_u <- nb_obs_g - K - 1

    if (nb_u >= 1) {
      DeltaY_X_g <- DeltaY_X[ind_g, , drop = FALSE]
      DeltaY_Xg_gamma <- DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 1] -
        DeltaY_X_g[(2 + K):nrow(DeltaY_X_g), 2:ncol(DeltaY_X_g), drop = FALSE] %*% gamma

      M_g <- M[(counter + 1):(counter + nb_u), , drop = FALSE]
      pinv_Mg <- MASS::ginv(t(M_g))

      for (param_idx in 1:n_params) {
        ek <- numeric(n_params)
        ek[param_idx] <- 1
        f_gk <- pinv_Mg %*% ek

        if (sqrt(sum((t(M_g) %*% f_gk - ek)^2)) < 1e-8) {
          ind_avg[j, param_idx] <- weights[j]
          beta_bar_hat[j, param_idx] <- as.numeric(t(f_gk) %*% DeltaY_Xg_gamma)
        }
      }
      counter <- counter + max(0, nb_obs_g - K - 1)
    }
    j <- j + 1
  }

  Nobs <- colSums(ind_avg)
  beta_bar_hat_final <- numeric(n_params)
  for (param_idx in 1:n_params) {
    if (Nobs[param_idx] > 0) {
      beta_bar_hat_final[param_idx] <- sum(ind_avg[, param_idx] * beta_bar_hat[, param_idx]) / Nobs[param_idx]
    }
  }

  # Create names for coefficients
  param_names <- character(n_params)
  idx <- 1
  # Main effects
  for (k in 0:K) {
    param_names[idx] <- paste0("beta_", k, ",", k)
    idx <- idx + 1
  }
  # Interaction effects
  for (k in 0:(K - 1)) {
    for (kp in (k + 1):K) {
      param_names[idx] <- paste0("beta_", k, ",", kp)
      idx <- idx + 1
    }
  }

  names(beta_bar_hat_final) <- param_names
  return(list(gamma = gamma, Nobs = Nobs, B_hat = beta_bar_hat_final, model = "interactions"))
}


#' Internal Bootstrap Function for RC Models
#'
#' Internal function to compute bootstrap standard errors and confidence intervals.
#' This function is called by estim_RC_model_unified when bootstrap = TRUE.
#'
#' @param result A result object from an RC model estimation function.
#' @param data Data frame containing the panel data.
#' @param group_col Character. Name of the column identifying groups.
#' @param deltaY_col Character. Name of the column for ΔY.
#' @param deltaD_col Character. Name of the column for ΔD.
#' @param D_col Character. Name of the column for D.
#' @param X_cols Character vector. Names of columns for covariates.
#' @param K Integer. Number of lags (required for "base" and "interactions" models).
#' @param B Integer. Number of bootstrap iterations.
#' @param conf_level Numeric. Confidence level for intervals.
#' @param verbose Logical. Show progress bar.
#'
#' @return A list containing bootstrap statistics:
#'   \item{se}{Bootstrap standard errors}
#'   \item{t_stat}{t-statistics}
#'   \item{ci_lower}{Lower bounds of confidence intervals}
#'   \item{ci_upper}{Upper bounds of confidence intervals}
#'   \item{boot_iterations}{Number of successful bootstrap iterations}
#'
#' @keywords internal
.bootstrap_RC_model_internal <- function(result, data, group_col, deltaY_col, deltaD_col,
                                         D_col, X_cols, K = NULL, B = 1000,
                                         conf_level = 0.95, verbose = TRUE) {

  # Validate inputs
  if (!is.list(result) || !all(c("B_hat", "model") %in% names(result))) {
    stop("result must be an output from an RC model estimation function")
  }

  model_type <- result$model

  # Validate K for models that need it
  if (model_type %in% c("base", "interactions") && is.null(K)) {
    stop("K must be specified for base and interactions models")
  }

  # Get unique groups
  list_groups <- unique(data[[group_col]])
  n_groups <- length(list_groups)

  # Storage for bootstrap results
  n_coef <- length(result$B_hat)
  Bhat_boot <- matrix(0, B, n_coef)

  # Progress indicator
  if (verbose) {
    cat(sprintf("Computing bootstrap standard errors (%d iterations)...\n", B))
    pb <- txtProgressBar(min = 0, max = B, style = 3)
  }

  # Bootstrap iterations
  for (b in 1:B) {
    # Sample groups with replacement
    groups_boot <- sample(list_groups, n_groups, replace = TRUE)

    # Calculate weights (how many times each group appears)
    weights <- as.numeric(table(factor(groups_boot, levels = list_groups)))

    # Estimate model on bootstrap sample
    tryCatch({
      if (model_type == "base") {
        result_boot <- estim_RC_model_base(
          K = K,
          data = data,
          group_col = group_col,
          deltaY_col = deltaY_col,
          deltaD_col = deltaD_col,
          D_col = D_col,
          X_cols = X_cols,
          weights = weights
        )
      } else if (model_type == "full_dynamics") {
        result_boot <- estim_RC_model_full(
          K = NULL,
          data = data,
          group_col = group_col,
          deltaY_col = deltaY_col,
          deltaD_col = deltaD_col,
          D_col = D_col,
          X_cols = X_cols,
          weights = weights
        )
      } else if (model_type == "interactions") {
        result_boot <- estim_RC_model_interactions(
          K = K,
          data = data,
          group_col = group_col,
          deltaY_col = deltaY_col,
          deltaD_col = deltaD_col,
          D_col = D_col,
          X_cols = X_cols,
          weights = weights
        )
      }

      Bhat_boot[b, ] <- result_boot$B_hat

    }, error = function(e) {
      Bhat_boot[b, ] <- NA
    })

    if (verbose) {
      setTxtProgressBar(pb, b)
    }
  }

  if (verbose) {
    close(pb)
  }

  # Remove failed iterations
  n_failed <- sum(!complete.cases(Bhat_boot))
  if (n_failed > 0) {
    warning(sprintf("%d bootstrap iterations failed and were removed", n_failed))
  }
  Bhat_boot <- Bhat_boot[complete.cases(Bhat_boot), , drop = FALSE]

  if (nrow(Bhat_boot) == 0) {
    stop("All bootstrap iterations failed")
  }

  # Compute standard errors
  se_boot <- apply(Bhat_boot, 2, sd)
  names(se_boot) <- names(result$B_hat)

  # Compute t-statistics
  t_stat <- as.numeric(result$B_hat) / se_boot
  names(t_stat) <- names(result$B_hat)

  # Compute confidence intervals
  alpha <- 1 - conf_level
  z_crit <- qnorm(1 - alpha / 2)
  ci_lower <- as.numeric(result$B_hat) - z_crit * se_boot
  ci_upper <- as.numeric(result$B_hat) + z_crit * se_boot
  names(ci_lower) <- names(result$B_hat)
  names(ci_upper) <- names(result$B_hat)

  if (verbose) {
    cat(sprintf("\nBootstrap completed: %d successful iterations (%.1f%% CI)\n",
                nrow(Bhat_boot), conf_level * 100))
  }

  # Return bootstrap statistics
  return(list(
    se = se_boot,
    t_stat = t_stat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    boot_iterations = nrow(Bhat_boot)
  ))
}
