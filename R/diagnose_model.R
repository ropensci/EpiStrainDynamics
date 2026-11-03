#' Diagnose model convergence and fit
#'
#' Checks MCMC convergence diagnostics --- [R-hat](https://mc-stan.org/docs/reference-manual/analysis.html)
#' and [effective sample size](https://mc-stan.org/docs/reference-manual/analysis.html) ---
#' against user-specified thresholds, and reports any parameters that fail
#' either check. A per-parameter table is also returned, labelling each
#' parameter with the pathogen and time point it refers to, so that
#' convergence problems can be located in time and by pathogen.
#'
#' @param fitted_model A fitted model object of class `EpiStrainDynamics.fit`,
#'   as returned by [fit_model()].
#' @param rhat_threshold R-hat threshold for convergence (default 1.1). Values
#'   above this threshold indicate chains have not mixed well.
#' @param eff_sample_threshold Effective sample size threshold (default 100).
#'   Values below this threshold indicate the posterior samples for a
#'   parameter are too autocorrelated to reliably estimate its distribution.
#'
#' @details
#' The `parameters` table identifies each parameter by the exact Stan name
#' before its brackets (`parameter_group`), and converts the bracketed
#' indices into a pathogen name and a time:
#'
#' * `log_incidence` is the time-varying log-incidence trajectory. Each
#'   element is matched to the pathogen and the observation time it belongs
#'   to.
#' * `b_spline_coef`, in penalised spline models, holds the spline
#'   coefficients. A coefficient does not belong to a single time point, so
#'   it is placed at its Greville abscissa (the mean of the interior knots its
#'   basis function spans) and converted to the time scale of the data.
#'   Coefficients near the edges of the spline can fall before the first or
#'   after the last observation. Some of these have no weight on any
#'   observation, so only the smoothing prior informs them; a flag on one of
#'   them does not mean the fitted trajectory (`log_incidence`) has failed to
#'   converge.
#' * `noisy_incidence` is the pathogen-specific incidence with noise, present
#'   when the model includes noise in individual pathogens. It is matched to
#'   pathogen and time in the same way as the trajectory.
#' * `tau` is matched to a pathogen when each pathogen has its own smoothing
#'   parameter (`"independent"` smoothing, or a single-pathogen model). With
#'   `"shared"` smoothing one `tau` applies to all pathogens, so `pathogen`
#'   is `NA`. With `"correlated"` smoothing there is no `tau`; the
#'   covariance matrix `Sigma` is used instead.
#' * All other parameters (e.g. `phi`, `eta`, `Sigma`,
#'   `day_of_week_simplex`, `lp__`) have `NA` for `pathogen` and `time`.
#'
#' @importFrom rstan summary
#'
#' @srrstats {BS4.3, BS5.3, BS5.4} implements a convergence checker
#'
#' @return An object of class `EpiStrainDynamics.diagnostics` (a list),
#'   returned invisibly, containing:
#'   \item{convergence}{Logical; `TRUE` if no parameter fails either
#'     threshold}
#'   \item{rhat_issues}{Character vector of parameter names with R-hat above
#'     `rhat_threshold`}
#'   \item{eff_sample_issues}{Character vector of parameter names with
#'     effective sample size below `eff_sample_threshold`}
#'   \item{max_rhat}{The largest R-hat value across all parameters}
#'   \item{min_neff}{The smallest effective sample size across all
#'     parameters}
#'   \item{summary}{The full `rstan::summary()` table the above are derived
#'     from}
#'   \item{parameters}{A data frame with one row per parameter and columns
#'     `parameter` (full Stan name), `parameter_group` (Stan name without
#'     indices), `pathogen`, `time_idx` (position on the observation-index
#'     scale, where 1 is the first observation; fractional for spline
#'     coefficients), `time` (on the time scale of the input data), `rhat`,
#'     `n_eff`, and the logical flags `rhat_issue` and `eff_sample_issue`.
#'     See Details.}
#'   \item{rhat_threshold, eff_sample_threshold}{The thresholds the checks
#'     used}
#' @export
#' @examplesIf interactive()
#' mod <- construct_model(
#'   pathogen_structure = single(
#'     data = sarscov2,
#'     case_timeseries = "cases",
#'     time = "date"
#'   ),
#'   method = random_walk()
#' )
#' fit <- fit_model(mod)
#' diag <- diagnose_model(fit)
#'
#' # Parameters that failed either check, with pathogen and time attached
#' diag$parameters[diag$parameters$rhat_issue | diag$parameters$eff_sample_issue, ]
diagnose_model <- function(fitted_model,
                           rhat_threshold = 1.1,
                           eff_sample_threshold = 100) {
  if (!inherits(fitted_model, "EpiStrainDynamics.fit")) {
    cli::cli_abort("{.arg fitted_model} must be an {.cls EpiStrainDynamics.fit} object, got {.cls {class(fitted_model)}}")
  }

  validate_rhat_threshold(rhat_threshold)
  validate_eff_sample_threshold(eff_sample_threshold)

  # Get model summary
  fit_summary <- rstan::summary(fitted_model$fit)$summary

  # Check R-hat values
  rhat_values <- fit_summary[, "Rhat"]
  rhat_issues <- names(rhat_values)[rhat_values > rhat_threshold & !is.na(rhat_values)]

  # Check effective sample sizes
  n_eff_values <- fit_summary[, "n_eff"]
  eff_sample_issues <- names(n_eff_values)[n_eff_values < eff_sample_threshold & !is.na(n_eff_values)]

  # Overall convergence assessment
  convergence <- length(rhat_issues) == 0 && length(eff_sample_issues) == 0

  # Additional diagnostics
  max_rhat <- max(rhat_values, na.rm = TRUE)
  min_neff <- min(n_eff_values, na.rm = TRUE)

  # Per-parameter diagnostics labelled by pathogen and time
  parameters <- label_parameter_diagnostics(
    fitted_model,
    fit_summary,
    rhat_threshold = rhat_threshold,
    eff_sample_threshold = eff_sample_threshold
  )

  diagnostics <- list(
    convergence = convergence,
    rhat_issues = rhat_issues,
    eff_sample_issues = eff_sample_issues,
    max_rhat = max_rhat,
    min_neff = min_neff,
    summary = fit_summary,
    parameters = parameters,
    rhat_threshold = rhat_threshold,
    eff_sample_threshold = eff_sample_threshold
  )
  class(diagnostics) <- c("EpiStrainDynamics.diagnostics", class(diagnostics))

  # Print summary
  cat("Model Convergence Diagnostics\n")
  cat("=============================\n")
  cat("Overall convergence:", ifelse(convergence, "GOOD", "ISSUES DETECTED"), "\n")
  cat("Maximum R-hat:", round(max_rhat, 3), "\n")
  cat("Minimum n_eff:", round(min_neff, 0), "\n")

  invisible(diagnostics)
}

#' Label Parameter Diagnostics by Pathogen and Time
#'
#' Builds a per-parameter table of R-hat and effective sample size, with the
#' bracketed Stan indices converted to pathogen names and times. Parameter
#' groups are matched on the full Stan name before the brackets, never on a
#' prefix.
#'
#' @param fitted_model A fitted model object of class `EpiStrainDynamics.fit`
#' @param fit_summary The `summary` matrix from `rstan::summary()`
#' @param rhat_threshold The R-hat convergence threshold
#' @param eff_sample_threshold The minimum effective sample size threshold
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#' @srrstats {G5.2a} every error statement is unique
#'
#' @return A data frame with one row per parameter
label_parameter_diagnostics <- function(fitted_model,
                                        fit_summary,
                                        rhat_threshold,
                                        eff_sample_threshold) {
  constructed_model <- fitted_model$constructed_model
  model_type <- class(constructed_model)[1]

  known_types <- c(
    "rw_single", "rw_multiple", "rw_subtyped",
    "ps_single", "ps_multiple", "ps_subtyped"
  )
  if (!model_type %in% known_types) {
    cli::cli_abort(c(
      "Cannot label parameters for model type {.val {model_type}}",
      "i" = "Expected one of {.val {known_types}}"
    ))
  }

  is_spline <- startsWith(model_type, "ps_")
  is_single <- endsWith(model_type, "_single")

  pathogen_names <- constructed_model$pathogen_names
  tsbl <- constructed_model$validated_tsbl
  obs_time <- tsbl[[tsibble::index_var(tsbl)]]
  num_data <- length(obs_time)

  par_names <- rownames(fit_summary)
  n_par <- length(par_names)
  par_group <- sub("\\[.*$", "", par_names)
  indices <- parse_parameter_indices(par_names)

  pathogen_idx <- rep(NA_integer_, n_par)
  time_idx <- rep(NA_real_, n_par)

  # Trajectories, spline coefficients and noise terms: [pathogen, index] in
  # multi-pathogen models, [index] in single-pathogen models. Rows with an
  # unexpected number of indices are left unlabelled.
  expected_n_idx <- if (is_single) 1L else 2L
  trajectory_rows <- par_group %in%
    c("log_incidence", "b_spline_coef", "noisy_incidence") &
    indices$n_idx == expected_n_idx

  if (is_single) {
    pathogen_idx[trajectory_rows] <- 1L
    time_idx[trajectory_rows] <- indices$first[trajectory_rows]
  } else {
    pathogen_idx[trajectory_rows] <- indices$first[trajectory_rows]
    time_idx[trajectory_rows] <- indices$second[trajectory_rows]
  }

  # Spline coefficients are positioned at their Greville abscissae
  if (is_spline) {
    coef_rows <- trajectory_rows & par_group == "b_spline_coef"
    greville <- greville_abscissae(
      constructed_model$standata$knots,
      constructed_model$standata$spline_degree
    )
    time_idx[coef_rows] <- greville[time_idx[coef_rows]]
  }

  # Smoothing parameters: one per pathogen only under "independent"
  # smoothing (cov_structure == 1); a single shared tau otherwise
  tau_rows <- par_group == "tau"
  if (is_single) {
    pathogen_idx[tau_rows] <- 1L
  } else if (identical(as.integer(constructed_model$standata$cov_structure), 1L)) {
    pathogen_idx[tau_rows] <- indices$first[tau_rows]
  }

  # Convert positions to the time scale of the data, keeping its class.
  # Whole-number positions within the data are looked up directly; others
  # (spline coefficients) are converted using the observation spacing.
  time <- obs_time[rep(NA_integer_, n_par)]
  exact_rows <- !is.na(time_idx) &
    time_idx == round(time_idx) &
    time_idx >= 1 &
    time_idx <= num_data
  time[exact_rows] <- obs_time[time_idx[exact_rows]]

  approx_rows <- !is.na(time_idx) & !exact_rows
  if (any(approx_rows)) {
    time[approx_rows] <- index_to_time(time_idx[approx_rows], obs_time)
  }

  rhat <- unname(fit_summary[, "Rhat"])
  n_eff <- unname(fit_summary[, "n_eff"])

  data.frame(
    parameter = par_names,
    parameter_group = par_group,
    pathogen = pathogen_names[pathogen_idx],
    time_idx = time_idx,
    time = time,
    rhat = rhat,
    n_eff = n_eff,
    rhat_issue = !is.na(rhat) & rhat > rhat_threshold,
    eff_sample_issue = !is.na(n_eff) & n_eff < eff_sample_threshold,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Parse Bracketed Indices from Stan Parameter Names
#'
#' Extracts up to two integer indices from names such as `log_incidence[2,15]`,
#' `tau[3]` or `phi`.
#'
#' @param par_names Character vector of Stan parameter names
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
#' @return A list with integer vectors `n_idx` (number of indices, 0 for
#'   unindexed names), `first` and `second` (`NA` where absent)
parse_parameter_indices <- function(par_names) {
  has_idx <- grepl("\\[[0-9,]+\\]$", par_names)

  inner <- rep(NA_character_, length(par_names))
  inner[has_idx] <- sub("^.*\\[([0-9,]+)\\]$", "\\1", par_names[has_idx])

  parts <- strsplit(inner, ",", fixed = TRUE)
  n_idx <- lengths(parts)
  n_idx[!has_idx] <- 0L

  list(
    n_idx = n_idx,
    first = vapply(parts, function(x) as.integer(x[1L]), integer(1)),
    second = vapply(
      parts,
      function(x) if (length(x) >= 2L) as.integer(x[2L]) else NA_integer_,
      integer(1)
    )
  )
}

#' Greville Abscissae of the Penalised Spline Basis
#'
#' Returns the position of each B-spline basis function on the
#' observation-index scale, as the mean of the interior knots it spans. The
#' extended knot sequence mirrors the `transformed data` block of the
#' `ps_*.stan` models: boundary knots repeated `spline_degree` times.
#'
#' @param knots Numeric vector of knots, as stored in `standata$knots`
#' @param spline_degree Polynomial degree of the spline
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
#' @return Numeric vector of length `length(knots) + spline_degree - 1`
greville_abscissae <- function(knots, spline_degree) {
  num_knots <- length(knots)
  ext_knots <- c(
    rep(knots[1L], spline_degree),
    knots,
    rep(knots[num_knots], spline_degree)
  )
  num_basis <- num_knots + spline_degree - 1L

  vapply(
    seq_len(num_basis),
    function(j) mean(ext_knots[(j + 1L):(j + spline_degree)]),
    numeric(1)
  )
}

#' Convert Observation-Index Positions to Data Time
#'
#' Maps positions on the observation-index scale (1 = first observation) to
#' the time scale of the data, assuming regular spacing as the models do.
#' Positions may be fractional or lie outside the observed range.
#'
#' Integer, `yearweek`, `yearmonth` and `yearquarter` times can only hold
#' whole units, so positions are rounded to the nearest unit (e.g. week);
#' `time_idx` keeps the exact position.
#'
#' @param idx Numeric vector of positions on the observation-index scale
#' @param obs_time The time index of the validated tsibble
#'
#' @noRd
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
#' @return A vector of times of the same class as `obs_time`, `NA` where
#'   conversion is not possible
index_to_time <- function(idx, obs_time) {
  missing_time <- obs_time[rep(NA_integer_, length(idx))]

  if (length(obs_time) < 2L) {
    return(missing_time)
  }

  conversion_failed <- function(cnd) {
    cli::cli_warn(c(
      "Could not convert spline coefficient positions to the time scale of the data",
      "x" = "{conditionMessage(cnd)}",
      "i" = "Use the {.field time_idx} column of {.field parameters} instead"
    ))
    missing_time
  }

  whole_units <- is.integer(obs_time) ||
    inherits(obs_time, c("yearweek", "yearmonth", "yearquarter"))

  tryCatch(
    {
      if (whole_units) {
        step <- as.numeric(obs_time[2L]) - as.numeric(obs_time[1L])
        time <- obs_time[1L] + as.integer(round((idx - 1) * step))
      } else {
        step <- obs_time[2L] - obs_time[1L]
        time <- obs_time[1L] + (idx - 1) * step
      }

      if (!identical(class(time), class(obs_time))) {
        cli::cli_abort(
          "Result has class {.cls {class(time)}}, not {.cls {class(obs_time)}}"
        )
      }
      time
    },
    error = conversion_failed,
    warning = conversion_failed
  )
}
