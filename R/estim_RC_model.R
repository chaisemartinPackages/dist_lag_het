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
#' @param B Integer. Number of bootstrap iterations (default: 200, ignored if bootstrap = FALSE).
#' @param conf_level Numeric. Confidence level for bootstrap intervals (default: 0.95, ignored if bootstrap = FALSE).
#' @param verbose Logical. Show bootstrap progress bar (default: TRUE, ignored if bootstrap = FALSE).
#' @param same_sample Logical. If TRUE, restrict estimation to groups identified for ALL
#'   coefficients simultaneously (default: FALSE). When TRUE, also returns \code{n_valid_groups}.
#'
#' @return A list containing:
#'   \item{gamma}{Coefficients for covariates}
#'   \item{B_hat}{Estimated treatment effect coefficients (betas)}
#'   \item{Nobs}{Number of observations used for each coefficient}
#'   \item{model}{Model type used}
#'   \item{n_valid_groups}{Number of groups identified for all coefficients (only if same_sample = TRUE)}
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
#'   B = 200
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
                                   bootstrap = FALSE, B = 200, conf_level = 0.95,
                                   verbose = TRUE, same_sample = FALSE) {

  model <- match.arg(model)

  # Estimate the main model
  if (model == "base") {
    if (is.null(K)) stop("K must be specified for base model")
    result <- estim_RC_model_base(K, data, group_col, deltaY_col, deltaD_col,
                               D_col, X_cols, weights, same_sample)

  } else if (model == "full_dynamics") {
    if (!is.null(K)) warning("K is ignored for full_dynamics model")
    result <- estim_RC_model_full(K = NULL, data, group_col, deltaY_col, deltaD_col,
                               D_col, X_cols, weights, same_sample)

  } else if (model == "interactions") {
    if (is.null(K)) stop("K must be specified for interactions model")
    result <- estim_RC_model_interactions(K, data, group_col, deltaY_col, deltaD_col,
                                       D_col, X_cols, weights, same_sample)
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
      verbose = verbose,
      same_sample = same_sample
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


#' Base RC Model Estimation (C++ backend)
#'
#' Estimates the base RC model with K lags using optimized C++ code.
#' This function uses RcppEigen for fast matrix operations.
#'
#' @inheritParams estim_RC_model_unified
#'
#' @return A list containing estimation results.
#'
#' @export
estim_RC_model_base <- function(K, data, group_col, deltaY_col, deltaD_col,
                                D_col, X_cols, weights = NULL, same_sample = FALSE) {

  # Convert to dataframe if needed
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Extract and convert to integer group IDs for faster C++ processing
  list_g_data <- as.integer(factor(data[[group_col]]))
  list_g <- sort(unique(list_g_data))
  nb_c <- length(list_g)

  # Extract numeric data
  DeltaY <- as.numeric(data[[deltaY_col]])
  DeltaD <- as.numeric(data[[deltaD_col]])
  D <- as.numeric(data[[D_col]])

  # Extract X as matrix
  if (length(X_cols) > 0) {
    X <- as.matrix(data[, X_cols, drop = FALSE])
    storage.mode(X) <- "double"
  } else {
    X <- matrix(0, nrow(data), 0)
  }

  # Set weights
  if (is.null(weights)) {
    weights <- rep(1, nb_c)
  }

  # Call C++ implementation
  result <- estim_RC_model_cpp_full(
    K = as.integer(K),
    list_g_data = list_g_data,
    list_g = list_g,
    DeltaY = DeltaY,
    DeltaD = DeltaD,
    D = D,
    X = X,
    weights = weights,
    same_sample = same_sample
  )

  # Add names and model type
  names(result$B_hat) <- paste0("beta_", 0:K)
  result$model <- "base"

  return(result)
}


#' Full Dynamics RC Model Estimation (C++ backend)
#'
#' Estimates the RC model with full dynamics (all available periods)
#' using optimized C++ code.
#'
#' @inheritParams estim_RC_model_unified
#'
#' @return A list containing estimation results.
#'
#' @export
estim_RC_model_full <- function(K = NULL, data, group_col, deltaY_col, deltaD_col,
                                D_col, X_cols, weights = NULL, same_sample = FALSE) {

  # Convert to dataframe if needed
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Extract and convert to integer group IDs
  list_g_data <- as.integer(factor(data[[group_col]]))
  list_g <- sort(unique(list_g_data))
  nb_c <- length(list_g)

  # Compute group sizes
  group_sizes <- as.integer(table(factor(list_g_data, levels = list_g)))

  # Determine T_periods
  T_periods <- max(group_sizes)

  # Extract numeric data
  DeltaY <- as.numeric(data[[deltaY_col]])
  DeltaD <- as.numeric(data[[deltaD_col]])

  # Extract X as matrix
  if (length(X_cols) > 0) {
    X <- as.matrix(data[, X_cols, drop = FALSE])
    storage.mode(X) <- "double"
  } else {
    X <- matrix(0, nrow(data), 0)
  }

  # Set weights
  if (is.null(weights)) {
    weights <- rep(1, nb_c)
  }

  # Call C++ implementation
  result <- estim_RC_model_full_cpp(
    T_periods = as.integer(T_periods),
    list_g_data = list_g_data,
    list_g = list_g,
    DeltaY = DeltaY,
    DeltaD = DeltaD,
    X = X,
    weights = weights,
    group_sizes = group_sizes,
    same_sample = same_sample
  )

  # Add names and model type
  names(result$B_hat) <- paste0("beta_", 0:(T_periods - 2))
  result$model <- "full_dynamics"

  return(result)
}


#' RC Model with Interaction Terms (C++ backend)
#'
#' Estimates the RC model with interaction terms between lags
#' using optimized C++ code.
#'
#' @inheritParams estim_RC_model_unified
#'
#' @return A list containing estimation results.
#'
#' @export
estim_RC_model_interactions <- function(K, data, group_col, deltaY_col, deltaD_col,
                                        D_col, X_cols, weights = NULL, same_sample = FALSE) {

  # Convert to dataframe if needed
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # Extract and convert to integer group IDs
  list_g_data <- as.integer(factor(data[[group_col]]))
  list_g <- sort(unique(list_g_data))
  nb_c <- length(list_g)

  # Extract numeric data
  DeltaY <- as.numeric(data[[deltaY_col]])
  DeltaD <- as.numeric(data[[deltaD_col]])
  D <- as.numeric(data[[D_col]])

  # Extract X as matrix
  if (length(X_cols) > 0) {
    X <- as.matrix(data[, X_cols, drop = FALSE])
    storage.mode(X) <- "double"
  } else {
    X <- matrix(0, nrow(data), 0)
  }

  # Set weights
  if (is.null(weights)) {
    weights <- rep(1, nb_c)
  }

  # Call C++ implementation
  result <- estim_RC_model_interactions_cpp(
    K = as.integer(K),
    list_g_data = list_g_data,
    list_g = list_g,
    DeltaY = DeltaY,
    DeltaD = DeltaD,
    D = D,
    X = X,
    weights = weights,
    same_sample = same_sample
  )

  # Create names for coefficients
  n_params <- (K + 1) * (K + 2) / 2
  param_names <- character(n_params)
  idx <- 1
  for (k in 0:K) {
    param_names[idx] <- paste0("beta_", k, ",", k)
    idx <- idx + 1
  }
  for (k in 0:(K - 1)) {
    for (kp in (k + 1):K) {
      param_names[idx] <- paste0("beta_", k, ",", kp)
      idx <- idx + 1
    }
  }

  names(result$B_hat) <- param_names
  result$model <- "interactions"

  return(result)
}


#' Internal Bootstrap Function for RC Models
#'
#' Internal function to compute bootstrap standard errors and confidence intervals.
#' This function is called by estim_RC_model_unified when bootstrap = TRUE.
#' Uses the fast C++ backend for base model.
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
                                         conf_level = 0.95, verbose = TRUE,
                                         same_sample = FALSE) {

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
        # Use C++ backend for speed
        result_boot <- estim_RC_model_base(
          K = K,
          data = data,
          group_col = group_col,
          deltaY_col = deltaY_col,
          deltaD_col = deltaD_col,
          D_col = D_col,
          X_cols = X_cols,
          weights = weights,
          same_sample = same_sample
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
          weights = weights,
          same_sample = same_sample
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
          weights = weights,
          same_sample = same_sample
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
