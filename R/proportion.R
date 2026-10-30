#' Generic Method for proportion Analysis
#'
#' Computes epidemiological proportion, defined as the relative fraction of
#' cases attributable to specific pathogen(s) or strain(s) at a given time point,
#' calculated as the ratio of incidence from selected pathogen(s) to a
#' reference group:
#' \deqn{P_t = \frac{I_{\text{numerator}}(t)}{I_{\text{denominator}}(t)}}
#'
#' Where incidences are derived from the exponential of log-incidence estimates:
#' \deqn{P_t = \frac{\sum_{i \in \text{numerator}} \exp(\log\text{-incidence}_{i,t})}{\sum_{j \in \text{denominator}} \exp(\log\text{-incidence}_{j,t})}}
#'
#' Where the numerator and denominator are each a user-specified set of
#' pathogens or subtypes (via `numerator_combination`/`denominator_combination`),
#' and \eqn{i}, \eqn{j} index the pathogens included in each.
#'
#' This metric quantifies the **relative contribution** of specific pathogen(s)
#' or strain(s) to the total disease burden. Key characteristics:
#' \itemize{
#'   \item Values between 0 and 1 (\eqn{0 \leq P_t \leq 1}) when denominator
#'    includes numerator components
#'   \item Values can exceed 1 (\eqn{P_t > 1}) when denominator excludes
#'   numerator components
#'   \item Represents the fractional share of cases at each time point
#'   \item Time-varying to capture changing pathogen/strain dynamics
#' }
#'
#' **Flexible combinations**:
#' \itemize{
#'   \item Individual proportions: Each pathogen relative to all pathogens
#'   (default)
#'   \item Custom numerator: Specific pathogen(s) of interest (e.g., variant
#'   of concern)
#'   \item Custom denominator: Either 'all' pathogens or a specified subset
#'   \item Subset comparisons: Compare specific groups (e.g., Alpha vs.
#'   Delta + Omicron)
#' }
#'
#' This metric function can be run directly on the fitted model output.
#'
#' @param fitted_model Fitted model object with class `EpiStrainDynamics.fit`
#'  with `multiple` or `subtyped` pathogen structure.
#' @param numerator_combination Named pathogens or subtypes to be included in
#'  proportion numerator, or NULL. If NULL, it will use each pathogen.
#' @param denominator_combination Named pathogens or subtypes to be included in
#'  proportion denominator, or NULL. If NULL, it will use all pathogens.
#' @param ... Additional arguments passed to metrics calculation
#' @return named list of class `EpiStrainDynamics.metric` containing a dataframe
#'  of the calculated metric outcome (`$measure`), the fit object (`$fit`), the
#'  constructed model object (`$constructed_model`), the resolved pathogen names
#'  used as the numerator (`$numerator_combination`), and the resolved pathogen
#'  names used as the denominator (`$denominator_combination`). When
#'  `numerator_combination` is left as the default (`NULL`), every pathogen name
#'  is listed in `$numerator_combination`, but each was computed as its own
#'  separate proportion line (one per pathogen), not summed into a single group
#'  the way `$denominator_combination` always is. The `measure` data frame
#'  contains the median of the epidemiological quantity (`y`), the 50% credible
#'  interval of the quantity (`lb_50` & `ub_50`), the 95% credible interval
#'  (`lb_95` & `ub_95`), the proportion greater than a defined threshold value
#'  (`prop`), the pathogen name (`pathogen`), and the time label (`time`).
#' @family metrics
#' @export
#'
#' @srrstats {G1.3} metric defined clearly
#' @srrstats {G1.4} uses `Roxygen2` documentation
#'
#' @examplesIf interactive()
#' mod <- construct_model(
#'   pathogen_structure = multiple(
#'     case_timeseries = sarscov2$cases,
#'     time = sarscov2$date,
#'     component_pathogen_timeseries = list(
#'       alpha = sarscov2$alpha,
#'       delta = sarscov2$delta,
#'       omicron = sarscov2$omicron,
#'       other = sarscov2$other
#'     )
#'   ),
#'   method = p_spline()
#' )
#'
#' fit <- fit_model(mod)
#' prop <- proportion(fit)
#'
#' # or a unique combination, compared to all pathogens
#' prop2 <- proportion(fit,
#'   numerator_combination = c("alpha", "delta", "omicron")
#' )
#'
#' # or a user-specified combination in both numerator and denominator
#' prop3 <- proportion(fit,
#'   numerator_combination = "alpha",
#'   denominator_combination = c("alpha", "delta", "omicron")
#' )
proportion <- function(fitted_model,
                       numerator_combination = NULL,
                       denominator_combination = NULL, ...) {
  validate_class_inherits(fitted_model, "EpiStrainDynamics.fit")
  validate_class_inherits(fitted_model, c("ps", "rw"), require_all = FALSE)

  # Validate combination arguments
  pathogen_names <- unique(fitted_model$constructed_model$pathogen_names)
  num_idx <- validate_pathogen_combination(
    numerator_combination, pathogen_names, "numerator_combination"
  )
  denom_idx <- validate_pathogen_combination(
    denominator_combination, pathogen_names, "denominator_combination"
  )

  use_splines <- ifelse(inherits(fitted_model, "ps"), TRUE, FALSE)

  if (is.null(numerator_combination)) {
    measure <- do.call(
      rbind, lapply(
        num_idx,
        function(x) {
          compute_proportion(
            fitted_model,
            use_splines = use_splines,
            numerator_idx = x,
            denominator_idx = denom_idx
          )
        }
      )
    )
  } else {
    measure <- compute_proportion(fitted_model,
      use_splines = use_splines,
      numerator_idx = num_idx,
      denominator_idx = denom_idx
    )
  }

  out <- list(
    measure = measure,
    fit = fitted_model$fit,
    constructed_model = fitted_model$constructed_model,
    numerator_combination = pathogen_names[num_idx],
    denominator_combination = pathogen_names[denom_idx]
  )
  class(out) <- c("proportion", "EpiStrainDynamics.metric", class(out))
  out
}

#' Generic Computation Engine for porportion calculation
#'
#' High-level function that coordinates analysis for single pathogen models
#'
#' @param fitted_model Fitted model object
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
compute_proportion <- function(fitted_model,
                               use_splines, ...) {
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

  extra_args <- list(...)

  # Create results
  selection_index <- 1:components$num_days
  time_grid <- data.frame(time_idx = selection_index)

  rep <- extra_args$numerator_idx # if not proportion calc then this will be NULL
  pathogen_idx_col <- rep(list(rep), nrow(time_grid))

  results <- calc_wrapper(time_grid, time_grid$time_idx,
    pathogen_idx_col = pathogen_idx_col,
    calc_proportion,
    a, post, components, extra_args, threshold = 0
  )

  measure <- cbind(results,
    time = components$time[selection_index]
  )


  measure$pathogen <- do.call(rbind, pathogen_idx_col)
  pathogen_names_vec <- components$pathogen_names
  measure$pathogen <- apply(
    measure$pathogen, 1,
    function(x) {
      paste(pathogen_names_vec[x],
        collapse = ", "
      )
    }
  )

  measure
}

#' Calculate proportion for Individual Pathogen
#'
#' Computes log-difference in incidence for a specific pathogen
#'
#' @param a Array of log-incidence posterior samples \code{[samples, pathogens, time]}
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#' @param numerator_idx Integer index of the pathogen/s used in numerator of proportion
#' @param denominator_idx Integer index of the pathogen/s used in denominator of proportion
#'
#' @return Vector of proportion posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_proportion <- function(a, time_idx, pathogen_idx, post, components,
                            numerator_idx, denominator_idx) {
  if (length(pathogen_idx) > 1) {
    num <- rowSums(exp(a[, pathogen_idx, time_idx]))
  } else {
    num <- exp(a[, pathogen_idx, time_idx])
  }

  den <- rowSums(exp(a[, denominator_idx, time_idx]))

  num / den
}
