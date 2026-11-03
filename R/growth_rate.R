#' Generic Method for Growth Rate Analysis
#'
#' Computes the epidemiological growth rate, defined as the instantaneous
#' rate of change in log-incidence over time. Mathematically, it represents:
#' \deqn{r_t = \log(I_t) - \log(I_{t-1}) = \log\left(\frac{I_t}{I_{t-1}}\right)}
#' where \eqn{I_t} is incidence (see [incidence()]) at time \eqn{t}.
#'
#' Growth rate and \eqn{R_t} both describe transmission trends, but on
#' different scales: growth rate is a direct log-scale rate of change, while
#' [Rt()] additionally accounts for the generation interval to translate
#' that rate into an average number of secondary infections per case.
#' Growth rate is positive whenever \eqn{R_t > 1} and negative whenever
#' \eqn{R_t < 1}, since both describe the same underlying growth or decline.
#'
#' This metric quantifies the **proportional change** in disease incidence
#' from one time period to the next on a logarithmic scale, where:
#' \itemize{
#'   \item Positive values (\eqn{r_t > 0}) indicate exponential growth
#'   \item Negative values (\eqn{r_t < 0}) indicate exponential decline
#'   \item Values near zero (\eqn{r_t \approx 0}) indicate stable incidence
#'   \item The magnitude indicates the rate of exponential change
#' }
#'
#' For example:
#' \itemize{
#'   \item A growth rate of 0.1 means incidence increased by approximately 10.5% (\eqn{e^{0.1} - 1})
#'   \item A growth rate of -0.05 means incidence decreased by approximately 4.9%
#'   \item A growth rate of 0 means no change in incidence
#' }
#'
#' This metric function can be run directly on the fitted model output.
#'
#' @param fitted_model Fitted model object with class `EpiStrainDynamics.fit`
#' @param ... Additional arguments passed to metrics calculation
#' @return named list of class `EpiStrainDynamics.metric` containing a dataframe
#'  of the calculated metric outcome (`$measure`), the fit object (`$fit`), and the
#'  constructed model object (`$constructed_model`). The `measure` data frame
#'  contains the median of the epidemiological quantity (`y`), the 50% credible
#'  interval of the quantity (`lb_50` & `ub_50`), the 95% credible interval
#'  (`lb_95` & `ub_95`), the posterior probability that the growth rate is
#'  above 0 (`prop`), the pathogen name (`pathogen`), and the time label
#'  (`time`).
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
#' gr <- growth_rate(fit)
growth_rate <- function(fitted_model, ...) {
  validate_class_inherits(fitted_model, "EpiStrainDynamics.fit")
  UseMethod("growth_rate")
}

#' @rdname growth_rate
#' @export
growth_rate.ps <- function(fitted_model, ...) {
  out <- compute_multi_pathogen(fitted_model, 2, "growth_rate",
    threshold = 0, use_splines = TRUE
  )
  class(out) <- c("growth_rate", "EpiStrainDynamics.metric", class(out))
  out
}

#' @rdname growth_rate
#' @export
growth_rate.rw <- function(fitted_model, ...) {
  out <- compute_multi_pathogen(fitted_model, 2, "growth_rate",
    threshold = 0, use_splines = FALSE
  )
  class(out) <- c("growth_rate", "EpiStrainDynamics.metric", class(out))
  out
}

#' @rdname growth_rate
#' @export
growth_rate.ps_single <- function(fitted_model, ...) {
  out <- compute_single_pathogen(fitted_model, 2, "growth_rate",
    threshold = 0, use_splines = TRUE
  )
  class(out) <- c("growth_rate", "EpiStrainDynamics.metric", class(out))
  out
}

#' @rdname growth_rate
#' @export
growth_rate.rw_single <- function(fitted_model, ...) {
  out <- compute_single_pathogen(fitted_model, 2, "growth_rate",
    threshold = 0, use_splines = FALSE
  )
  class(out) <- c("growth_rate", "EpiStrainDynamics.metric", class(out))
  out
}

# =====================
# CALCULATION FUNCTIONS
# =====================

#' Calculate Growth Rate for Single Pathogen Model
#'
#' Computes log-difference in incidence for single pathogen models
#'
#' @param a Array of log-incidence posterior samples \code{[samples, time]}
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object (unused but required for interface
#'   consistency)
#' @param components Model components (unused but required for interface
#'   consistency)
#'
#' @return Vector of growth rate posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_growth_single <- function(a, time_idx, pathogen_idx, post, components) {
  a[, time_idx] - a[, time_idx - 1]
}

#' Calculate Growth Rate for Individual Pathogen
#'
#' Computes log-difference in incidence for a specific pathogen
#'
#' @param a Array of log-incidence posterior samples \code{[samples, pathogens,
#'  time]}
#' @param time_idx Integer time index
#' @param pathogen_idx Integer pathogen index
#' @param post Posterior samples object (unused but required for interface
#'  consistency)
#' @param components Model components (unused but required for interface
#' consistency)
#'
#' @return Vector of growth rate posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_growth_individual <- function(a, time_idx, pathogen_idx, post, components) {
  a[, pathogen_idx, time_idx] - a[, pathogen_idx, time_idx - 1]
}

#' Calculate Total Growth Rate Across All Pathogens
#'
#' Computes log-difference in total incidence across all pathogens
#'
#' @param a Array of log-incidence posterior samples \code{[samples, pathogens,
#'  time]}
#' @param time_idx Integer time index
#' @param pathogen_idx NULL (unused but required for interface consistency)
#' @param post Posterior samples object (unused but required for interface consistency)
#' @param components Model components (unused but required for interface consistency)
#'
#' @return Vector of total growth rate posterior samples
#' @noRd
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.4a} internal function specified with `@noRd`
#'
calc_growth_total <- function(a, time_idx, pathogen_idx, post, components) {
  total_i <- rowSums(exp(a[, , time_idx]))
  total_im1 <- rowSums(exp(a[, , time_idx - 1]))
  log(total_i) - log(total_im1)
}
