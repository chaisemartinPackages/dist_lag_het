#' Generate Synthetic DID Panel Data
#'
#' Creates synthetic panel data for testing and examples. The data includes
#' groups, time periods, treatment indicators, outcomes, and covariates.
#'
#' @param n_groups Integer. Number of groups (e.g., counties, individuals). Default is 50.
#' @param n_periods Integer. Number of time periods. Default is 5.
#' @param treatment_prob Numeric. Probability of treatment switching (0 to 1). Default is 0.3.
#' @param seed Integer. Random seed for reproducibility. Default is 123.
#'
#' @return A data frame with columns:
#'   \item{group}{Group identifier}
#'   \item{time}{Time period}
#'   \item{D}{Treatment level}
#'   \item{delta_D}{First difference of treatment}
#'   \item{Y}{Outcome variable}
#'   \item{delta_Y}{First difference of outcome}
#'   \item{X1, X2}{Covariates}
#'
#' @examples
#' # Generate data with 30 groups and 4 time periods
#' data <- generate_did_data(n_groups = 30, n_periods = 4)
#' head(data)
#'
#' @export
generate_did_data <- function(n_groups = 50, n_periods = 5, treatment_prob = 0.3, seed = 123) {
  set.seed(seed)

  # Initialize storage
  data_list <- list()

  for (g in 1:n_groups) {
    # Generate treatment pattern
    D <- numeric(n_periods)

    # Randomly decide if this group will ever be treated
    if (runif(1) < treatment_prob) {
      # Choose a random period to start treatment (not in first period)
      treatment_start <- sample(2:n_periods, 1)
      D[treatment_start:n_periods] <- 1
    }

    # Generate covariates (constant within group over time for simplicity)
    X1 <- rnorm(1, mean = 0, sd = 1)
    X2 <- rnorm(1, mean = 0, sd = 1)

    # Generate outcome
    # Y depends on treatment, covariates, group fixed effect, and time fixed effect
    group_effect <- rnorm(1, mean = 0, sd = 2)

    Y <- numeric(n_periods)
    for (t in 1:n_periods) {
      time_effect <- 0.5 * t  # Linear time trend
      treatment_effect <- 2 * D[t]  # Treatment effect

      # Add dynamic treatment effects (if treated in previous period)
      if (t > 1 && D[t - 1] == 1) {
        treatment_effect <- treatment_effect + 1.5  # Lagged effect
      }

      Y[t] <- group_effect + time_effect + treatment_effect +
              0.8 * X1 + 1.2 * X2 + rnorm(1, mean = 0, sd = 1)
    }

    # Calculate first differences
    delta_D <- c(NA, diff(D))
    delta_Y <- c(NA, diff(Y))

    # Create group data frame
    group_data <- data.frame(
      group = g,
      time = 1:n_periods,
      D = D,
      delta_D = delta_D,
      Y = Y,
      delta_Y = delta_Y,
      X1 = X1,
      X2 = X2
    )

    data_list[[g]] <- group_data
  }

  # Combine all groups
  data <- do.call(rbind, data_list)
  rownames(data) <- NULL

  # Remove first period (NAs in differences)
  data <- data[!is.na(data$delta_Y), ]

  return(data)
}


#' Print Summary of RC Model Results
#'
#' Provides a formatted summary of estimation results, including bootstrap
#' standard errors and confidence intervals if available.
#'
#' @param result A list returned by estim_RC_model_unified or related functions.
#'
#' @return Invisibly returns the result object.
#'
#' @examples
#' \dontrun{
#' data <- generate_did_data()
#' result <- estim_RC_model_unified(K = 2, data = data, group_col = "group",
#'                                   deltaY_col = "delta_Y", deltaD_col = "delta_D",
#'                                   D_col = "D", X_cols = c("X1", "X2"),
#'                                   model = "base", bootstrap = TRUE)
#' summary_RC_model(result)
#' }
#'
#' @export
summary_RC_model <- function(result) {
  cat("\n========================================\n")
  cat("RC Model Estimation Results\n")
  cat("========================================\n\n")

  cat("Model type:", result$model, "\n")

  # Check if bootstrap was performed
  has_bootstrap <- !is.null(result$se)

  if (has_bootstrap) {
    cat(sprintf("Bootstrap iterations: %d\n", result$boot_iterations))
  }
  cat("\n")

  cat("Covariate coefficients (gamma):\n")
  print(result$gamma)
  cat("\n")

  cat("Treatment effect coefficients (B_hat):\n")

  # Create data frame with appropriate columns
  if (has_bootstrap) {
    coef_df <- data.frame(
      Coefficient = names(result$B_hat),
      Estimate = as.numeric(result$B_hat),
      Std_Error = as.numeric(result$se),
      t_stat = as.numeric(result$t_stat),
      CI_lower = as.numeric(result$ci_lower),
      CI_upper = as.numeric(result$ci_upper),
      N_obs = result$Nobs
    )
  } else {
    coef_df <- data.frame(
      Coefficient = names(result$B_hat),
      Estimate = as.numeric(result$B_hat),
      N_obs = result$Nobs
    )
  }

  print(coef_df, row.names = FALSE)
  cat("\n")

  if (has_bootstrap) {
    cat("Note: Standard errors and confidence intervals computed via bootstrap.\n\n")
  }

  invisible(result)
}
