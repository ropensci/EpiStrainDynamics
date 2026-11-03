#' Generic Method for Incidence Analysis
#'
#' Computes epidemiological incidence, defined as the number of new cases
#' occurring at a specific time point, derived by exponentiating the
#' log-incidence estimates from the fitted model:
#' \deqn{I_t = \exp(\log\text{-incidence}_t)}
#'
#' This metric quantifies the **absolute number of new cases** at each time
#' point, where it is:
#' \itemize{
#'   \item Always positive (since it's an exponentiated value)
#'   \item Represents the expected case count at time t
#'   \item Can be adjusted for day-of-week effects when modeled
#'   \item Provides uncertainty quantification through posterior credible intervals
#' }
#'
#' **Day-of-week adjustment**: When day-of-week effects are included in the model,
#' the incidence is further adjusted as:
#' \deqn{I_t^{adj} = I_t \times \text{week\_effect} \times \text{dow\_simplex}[\text{DOW}(t)]}
#'
#' Where:
#' \itemize{
#'   \item \eqn{\text{week\_effect}} is the number of distinct days modelled
#'   (7 for a full weekly cycle)
#'   \item \eqn{\text{dow\_simplex}} gives the relative reporting weight for
#'   each day of the week, estimated from the data
#'   \item \eqn{\text{DOW}(t)} maps time \eqn{t} to its day of the week
#' }
#' This accounts for systematic variations in case reporting (e.g. lower
#' weekend reporting) that are not part of the underlying transmission trend.
#'
#' This metric function can be run directly on the fitted model output.
#'
#' @param fitted_model Fitted model object with class `EpiStrainDynamics.fit`
#' @param dow Logical indicating whether to include day-of-week effects.
#'   If \code{NULL} or \code{NA} (default), uses the day-of-week setting from
#'   the fitted model. If \code{TRUE}, includes day-of-week effects (model must
#'   have been fitted with \code{dow_effect = TRUE}). If \code{FALSE}, excludes
#'   day-of-week effects.
#' @param ... Additional arguments passed to metrics calculation
#'
#' @return named list of class `EpiStrainDynamics.metric` containing a dataframe
#'  of the calculated metric outcome (`$measure`), the fit object (`$fit`), and the
#'  constructed model object (`$constructed_model`). The `measure` data frame
#'  contains the median of the epidemiological quantity (`y`), the 50% credible
#'  interval of the quantity (`lb_50` & `ub_50`), the 95% credible interval
#'  (`lb_95` & `ub_95`), the pathogen name (`pathogen`), and the time label
#'  (`time`).
#'
#' @family metrics
#' @export
#'
#' @srrstats {G1.3} metric defined clearly
#' @srrstats {G1.4} uses `Roxygen2` documentation
#'
#' @examplesIf interactive()
#' mod <- construct_model(
#'   pathogen_structure = single(
#'     case_timeseries = sarscov2$cases,
#'     time = sarscov2$date
#'   ),
#'   method = random_walk()
#' )
#'
#' fit <- fit_model(mod)
#'
#' # Use model's dow setting (default)
#' inc <- incidence(fit)
#'
#' # Explicitly exclude dow effects
#' inc_no_dow <- incidence(fit, dow = FALSE)
#'
#' # Explicitly include dow effects (if model has them)
#' inc_with_dow <- incidence(fit, dow = TRUE)
incidence <- function(fitted_model, dow = NULL, ...) {
  validate_class_inherits(fitted_model, "EpiStrainDynamics.fit")
  UseMethod("incidence")
}

#' @rdname incidence
#' @export
incidence.ps <- function(fitted_model, dow = NULL, ...) {
  # Resolve and validate dow
  dow <- resolve_dow(fitted_model, dow)

  out <- compute_multi_pathogen(fitted_model, 1, "incidence",
    threshold = NULL, use_splines = TRUE, dow
  )
  class(out) <- c("incidence", "EpiStrainDynamics.metric", class(out))
  out
}

#' @rdname incidence
#' @export
incidence.rw <- function(fitted_model, dow = NULL, ...) {
  # Resolve and validate dow
  dow <- resolve_dow(fitted_model, dow)

  out <- compute_multi_pathogen(fitted_model, 1, "incidence",
    threshold = NULL, use_splines = FALSE, dow
  )
  class(out) <- c("incidence", "EpiStrainDynamics.metric", class(out))
  out
}

#' @rdname incidence
#' @export
incidence.ps_single <- function(fitted_model, dow = NULL, ...) {
  # Resolve and validate dow
  dow <- resolve_dow(fitted_model, dow)

  out <- compute_single_pathogen(fitted_model, 1, "incidence",
    threshold = NULL, use_splines = TRUE, dow
  )
  class(out) <- c("incidence", "EpiStrainDynamics.metric", class(out))
  out
}

#' @rdname incidence
#' @export
incidence.rw_single <- function(fitted_model, dow = NULL, ...) {
  # Resolve and validate dow
  dow <- resolve_dow(fitted_model, dow)

  out <- compute_single_pathogen(fitted_model, 1, "incidence",
    threshold = NULL, use_splines = FALSE, dow
  )
  class(out) <- c("incidence", "EpiStrainDynamics.metric", class(out))
  out
}

#' Resolve and validate dow argument
#'
#' Helper function to handle dow argument defaults and validation.
#' Used by metric functions to ensure consistent dow handling.
#'
#' @param fitted_model Fitted model object with dow_effect in constructed_model
#' @param dow User-provided dow value (NULL, NA, TRUE, or FALSE)
#'
#' @return Logical value for dow (TRUE or FALSE)
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
resolve_dow <- function(fitted_model, dow) {
  # Validate type
  if (!is.null(dow) && !is.logical(dow)) {
    cli::cli_abort(
      "{.arg dow} must be {.code NULL} or a logical value ({.code TRUE}/{.code FALSE})."
    )
  }

  # Set default if NULL or NA
  if (is.null(dow) || isTRUE(is.na(dow))) {
    dow <- fitted_model$constructed_model$dow_effect
    cli::cli_inform(
      "Using day-of-week setting from model: {.val {dow}}"
    )
  }

  # Validate compatibility
  if (dow && !fitted_model$constructed_model$dow_effect) {
    cli::cli_abort(
      "Day-of-week effects cannot be incorporated because the model was not
      fitted with {.code dow_effect = TRUE}."
    )
  }

  dow
}

#' Calculate Incidence for Single Pathogen Model (with DOW adjustment)
#'
#' Computes incidence for single pathogen models, applying day-of-week effects if present
#'
#' @param a Array of log-incidence posterior samples \code{[samples, time]}
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object containing day_of_week_simplex if DOW effects used
#' @param components Model components containing DOW information
#' @param dow_effect Logical whether day-of-week effect to be calculated
#'
#' @return Vector of incidence posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_incidence_single <- function(a, time_idx, pathogen_idx,
                                  post, components, dow_effect) {
  incidence <- exp(a[, time_idx])

  # Apply DOW effect if required
  if (dow_effect) {
    incidence <- incidence * components$week_effect *
      post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}

#' Calculate Incidence for Individual Pathogen (with DOW adjustment)
#'
#' Computes incidence for a specific pathogen, applying day-of-week effects if present
#'
#' @param a Array of log-incidence posterior samples \code{[samples, pathogens, time]}
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object containing day_of_week_simplex if DOW effects used
#' @param components Model components containing DOW information
#' @param dow_effect Logical whether day-of-week effect to be calculated
#'
#' @return Vector of incidence posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_incidence_individual <- function(a, time_idx, pathogen_idx,
                                      post, components, dow_effect) {
  incidence <- exp(a[, pathogen_idx, time_idx])

  # Apply DOW effect if required
  if (dow_effect) {
    incidence <- incidence * components$week_effect *
      post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}

#' Calculate Total Incidence Across All Pathogens (with DOW adjustment)
#'
#' Computes total incidence across all pathogens, applying day-of-week effects if present
#'
#' @param a Array of log-incidence posterior samples \code{[samples, pathogens, time]}
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object containing day_of_week_simplex if DOW effects used
#' @param components Model components containing DOW information
#' @param dow_effect Logical whether day-of-week effect to be calculated
#'
#' @return Vector of total incidence posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_incidence_total <- function(a, time_idx, pathogen_idx,
                                 post, components, dow_effect) {
  incidence <- rowSums(exp(a[, , time_idx]))

  # Apply DOW effect if required
  if (dow_effect) {
    incidence <- incidence * components$week_effect *
      post$day_of_week_simplex[, components$DOW[time_idx]]
  }

  incidence
}
