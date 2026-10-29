#' Calculate Summary Statistics for Posterior Samples
#'
#' Computes quantiles and proportion above threshold for a vector of posterior samples
#'
#' @param values Numeric vector of posterior samples
#' @param threshold Numeric threshold for calculating proportion above (default: 0)
#'
#' @importFrom stats quantile
#'
#' @return Data frame with columns: y (median), lb_50/ub_50 (50% CI),
#'   lb_95/ub_95 (95% CI), prop (proportion > threshold)
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_stats <- function(values, threshold = 0) {
  quan <- stats::quantile(values, c(0.5, 0.025, 0.25, 0.75, 0.975))
  prop <- mean(values > threshold)

  data.frame(
    y = quan[[1]],
    lb_50 = quan[[3]],
    ub_50 = quan[[4]],
    lb_95 = quan[[2]],
    ub_95 = quan[[5]],
    prop = prop
  )
}

#' Extract Model Components from Fitted Model
#'
#' Extracts commonly used components from a fitted model object for analysis functions
#'
#' @param fitted_model Fitted model object containing fit, constructed_model, etc.
#'
#' @return List containing extracted components: fit, pathogen_names, num_path,
#'   time_seq, time, num_days, days_per_knot, spline_degree, DOW, week_effect,
#'   dow_effect
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
get_model_components <- function(fitted_model) {
  tsbl <- fitted_model$constructed_model$validated_tsbl
  list(
    fit = fitted_model$fit,
    pathogen_names = fitted_model$constructed_model$pathogen_names %||% NULL,
    num_path = length(fitted_model$constructed_model$pathogen_names %||% 1),
    time = tsbl[[tsibble::index(tsbl)]],
    num_days = nrow(tsbl),
    knots = fitted_model$constructed_model$standata$knots %||% NULL,
    spline_degree = fitted_model$constructed_model$standata$spline_degree %||% NULL,
    DOW = fitted_model$constructed_model$standata$DOW %||% NULL,
    week_effect = fitted_model$constructed_model$standata$week_effect %||% NULL,
    dow_effect = fitted_model$constructed_model$dow_effect
  )
}

#' Expand Time Grid for Multiple Pathogens
#'
#' Expands a time grid to include all pathogen-time combinations (sequential
#' pairing, not Cartesian product)
#'
#' @param time_grid Data frame from with time_idx columns
#' @param pathogen_names Character vector of pathogen names
#'
#' @return Data frame with columns: pathogen, pathogen_idx, time_idx
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
expand_pathogen_grid <- function(time_grid, pathogen_names) {
  num_days <- nrow(time_grid)
  num_path <- length(pathogen_names)

  data.frame(
    pathogen = rep(pathogen_names, each = num_days),
    pathogen_idx = rep(seq_along(pathogen_names), each = num_days),
    time_idx = rep(time_grid$time_idx, times = num_path)
  )
}

#' Predict B-spline Basis Matrix
#'
#' Creates B-spline basis matrix for transforming spline coefficients to time series
#'
#' @param time_seq Numeric vector of time points
#' @param days_per_knot Integer number of days between knots
#' @param spline_degree Integer degree of B-splines
#'
#' @keywords internal
#'
#' @importFrom splines bs
#' @importFrom stats predict
#'
#' @return Matrix B_true for transforming spline coefficients
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
predict_B_true <- function(time_seq, knots, spline_degree) {
  X <- as.numeric(time_seq)
  B_true <- splines::bs(seq(knots[1], knots[length(knots)], 1),
    knots = knots[2:(length(knots) - 1)],
    degree = spline_degree,
    intercept = TRUE
  )
  t(stats::predict(B_true, X))
}

#' Transform Posterior Samples Using Splines (Single Pathogen)
#'
#' Transforms spline coefficients to time series for single pathogen models
#'
#' @param post Posterior samples object containing spline coefficients
#' @param B_true B-spline basis matrix from predict_B_true()
#' @param num_days Integer number of time points
#'
#' @return Array of transformed posterior samples \code{[samples, time]}
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
transform_posterior_single <- function(post, B_true, num_days) {
  a <- array(data = NA, dim = c(nrow(post$b_spline_coef), num_days))
  for (k in seq_len(nrow(post$b_spline_coef))) {
    a[k, ] <- as.array(post$b_spline_coef[k, ]) %*% B_true
  }
  a
}

#' Transform Posterior Samples Using Splines (Multiple Pathogens)
#'
#' Transforms spline coefficients to time series for multiple pathogen models
#'
#' @param post Posterior samples object containing spline coefficients
#' @param B_true B-spline basis matrix from predict_B_true()
#' @param num_path Integer number of pathogens
#' @param num_days Integer number of time points
#'
#' @return Array of transformed posterior samples \code{[samples, pathogens, time]}
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
transform_posterior_multi <- function(post, B_true, num_path, num_days) {
  a <- array(data = NA, dim = c(nrow(post$b_spline_coef), num_path, num_days))
  for (j in 1:num_path) {
    for (k in seq_len(nrow(post$b_spline_coef))) {
      a[k, j, ] <- as.array(post$b_spline_coef[k, j, ]) %*% B_true
    }
  }
  a
}

#' Generic Computation Engine for Multiple Pathogen Analysis
#'
#' High-level function that coordinates analysis for multiple pathogen models
#'
#' @param fitted_model Fitted model object
#' @param start_idx Integer starting time index for analysis
#' @param measure Character string specifying metric ("incidence", "growth_rate", "Rt")
#' @param threshold Numeric threshold for proportion calculations (default: 0)
#' @param use_splines Logical indicating whether to use spline transformation
#' @param ... Additional arguments passed to calculation functions
#'  (e.g., tau_max, gi_dist for Rt, and dow for incidence)
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
#' @return Data frame with results for individual pathogens and totals
#' @importFrom rstan extract
#' @importFrom dplyr bind_rows arrange
#'
compute_multi_pathogen <- function(fitted_model, start_idx, measure,
                                   threshold = 0, use_splines = FALSE,
                                   ...) {
  components <- get_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    time_seq <- seq_len(components$num_days)
    B_true <- predict_B_true(
      time_seq, components$knots,
      components$spline_degree
    )
    a <- transform_posterior_multi(
      post, B_true, components$num_path,
      components$num_days
    )
  } else {
    a <- post$log_incidence
  }

  selection_index <- start_idx:components$num_days
  time_grid <- data.frame(time_idx = selection_index)
  extra_args <- list(...)

  # Individual pathogen results
  pathogen_grid <- expand_pathogen_grid(time_grid, components$pathogen_names)

  calc_individual_pathogen_fn <- switch(measure,
    "incidence" = calc_incidence_individual,
    "growth_rate" = calc_growth_individual,
    "Rt" = calc_rt_individual,
    {
      cli::cli_abort("Invalid option provided: '{measure}'.")
    }
  )
  pathogen_results <- calc_wrapper(
    pathogen_grid, pathogen_grid$time_idx,
    pathogen_grid$pathogen_idx,
    calc_individual_pathogen_fn,
    a, post, components, extra_args, threshold
  )

  calc_total_pathogens_fn <- switch(measure,
    "incidence" = calc_incidence_total,
    "growth_rate" = calc_growth_total,
    "Rt" = calc_rt_total,
    {
      cli::cli_abort("Invalid option provided: '{measure}'.")
    }
  )
  total_results <- calc_wrapper(time_grid, time_grid$time_idx,
    pathogen_idx_col = rep(list(NULL), nrow(time_grid)),
    calc_total_pathogens_fn,
    a, post, components, extra_args, threshold
  )
  total_results$pathogen <- "Total"

  measure_out <- dplyr::bind_rows(pathogen_results, total_results) |>
    dplyr::arrange(.data$pathogen != "Total", .data$pathogen) |>
    cbind(time = components$time[selection_index])

  #' @srrstats {TS4.0, TS4.0b, TS4.2, TS4.3} Return values are in unified class,
  #' which is documented.
  tsbl_measure <- tsibble::as_tsibble(
    measure_out,
    index = "time",
    key = "pathogen"
  )

  list(
    measure = tsbl_measure,
    fit = fitted_model$fit,
    constructed_model = fitted_model$constructed_model
  )
}

#' Generic Computation Engine for Single Pathogen Analysis
#'
#' High-level function that coordinates analysis for single pathogen models
#'
#' @param fitted_model Fitted model object
#' @param start_idx Integer starting time index for analysis
#' @param measure Character string specifying metric ("incidence", "growth_rate", "Rt")
#' @param threshold Numeric threshold for proportion calculations (default: 0)
#' @param use_splines Logical indicating whether to use spline transformation
#' @param ... Additional arguments passed to calculation functions (e.g., tau_max, gi_dist for Rt)
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
#' @return Data frame with analysis results
#' @importFrom rstan extract
#'
compute_single_pathogen <- function(fitted_model, start_idx, measure,
                                    threshold = 0, use_splines = FALSE, ...) {
  components <- get_model_components(fitted_model)
  post <- rstan::extract(components$fit)

  # Transform data if using splines
  if (use_splines) {
    time_seq <- seq_len(components$num_days)
    B_true <- predict_B_true(time_seq, components$knots, components$spline_degree)
    a <- transform_posterior_single(post, B_true, components$num_days)
  } else {
    a <- post$log_incidence
  }

  extra_args <- list(...)

  # Create results
  selection_index <- start_idx:components$num_days
  time_grid <- data.frame(time_idx = selection_index)

  calc_single_pathogen_fn <- switch(measure,
    "incidence" = calc_incidence_single,
    "growth_rate" = calc_growth_single,
    "Rt" = calc_rt_single,
    {
      cli::cli_abort("Invalid option provided: '{measure}'.")
    }
  )

  pathogen_idx_col <- rep(list(NULL), nrow(time_grid))

  results <- calc_wrapper(time_grid, time_grid$time_idx,
    pathogen_idx_col = pathogen_idx_col,
    calc_single_pathogen_fn,
    a, post, components, extra_args, threshold
  )

  results$pathogen <- components$pathogen_names

  #' @srrstats {TS4.0, TS4.0b, TS4.2, TS4.3} Return values are in unified class,
  #' which is documented.
  tsbl_measure <- tsibble::tsibble(
    time = components$time[selection_index],
    results,
    index = "time"
  )

  out <- list(
    measure = tsbl_measure,
    fit = fitted_model$fit,
    constructed_model = fitted_model$constructed_model
  )
  out
}

#' Unified Calculation Wrapper
#'
#' Wrapper function that applies calculation functions to data and computes summary statistics
#'
#' @param df Data frame containing time grid or pathogen grid
#' @param time_idx_col Vector of time indices
#' @param pathogen_idx_col Vector of pathogen indices (or list of NULLs)
#' @param calc_fn Calculation function to apply
#' @param a Array of posterior samples
#' @param post Posterior samples object
#' @param components Model components from get_model_components()
#' @param extra_args List of additional arguments for calc_fn
#' @param threshold Numeric threshold for summary statistics
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
#' @return Data frame with expanded summary statistics
#' @importFrom purrr map2
calc_wrapper <- function(df, time_idx_col, pathogen_idx_col, calc_fn,
                         a, post, components, extra_args, threshold) {
  # Calculate statistics for each time/pathogen combination
  stats_list <- purrr::map2(
    time_idx_col, pathogen_idx_col,
    function(t_idx, p_idx) {
      values <- do.call(
        calc_fn, c(
          list(
            a = a,
            time_idx = t_idx,
            pathogen_idx = p_idx,
            post = post,
            components = components
          ),
          extra_args
        )
      )
      calc_stats(values, threshold = threshold)
    }
  )

  # Convert list of dataframes to single data frame
  stats_df <- do.call(rbind, stats_list)

  # Combine with original data frame (excluding time_idx column)
  result_df <- cbind(df[, !names(df) %in% "time_idx", drop = FALSE], stats_df)

  result_df
}
