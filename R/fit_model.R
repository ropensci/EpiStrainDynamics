#' Generic Method for fitting model
#'
#' S3 generic for fitted models from constructed model object
#'
#' @param constructed_model prepared model object of class
#'  `EpiStrainDynamics.model`
#' @param n_chain number of MCMC chains, defaults to 4
#' @param n_iter A positive integer specifying the number of iterations for each
#'   chain, default value is 2000
#' @param n_warmup A positive integer specifying the number of warmup
#'   iterations,default value is half the number of iterations
#' @param thin A positive integer specifying the period for saving samples,
#'   default value is 1.
#' @param adapt_delta Numeric value between 0 and 1 indicating target average
#'   acceptance probability used in `rstan::sampling`. Default value is 0.9.
#' @param multi_cores A logical value indicating whether to parallelize chains
#'   with multiple cores, default is TRUE and uses all available cores - 1.
#' @param verbose Logical value controlling the verbosity of output. When TRUE
#'   (default), shows all messages, warnings, errors, and progress indicators.
#'   When FALSE, suppresses messages and progress while retaining warnings and
#'   errors.
#' @param suppress_warnings Logical value indicating whether to suppress
#'   warnings from Stan. Default is FALSE. When TRUE, warnings are suppressed
#'   but errors are still raised.
#' @param seed A positive integer seed used for random number generation in
#'   MCMC. Default is NULL, which means the seed is generated from 1 to the
#'   maximum integer supported by R.
#' @param ... additional arguments to `rstan::sampling()`, such as `init`
#'
#' @returns fit model of class `EpiStrainDynamics.fit`, or if fitting fails,
#'   an error is raised that can be caught and inspected. Its `fit` element is
#'   the `stanfit` object, whose main parameters are:
#'   * `log_incidence`: the latent log-incidence trend, indexed `[time]` for
#'     single pathogen models and `[pathogen, time]` otherwise. In p-spline
#'     models it is derived from `b_spline_coef`.
#'   * `b_spline_coef`: the B-spline coefficients (p-spline models only).
#'   * `noisy_incidence`: per-pathogen incidence with gamma noise around
#'     `exp(log_incidence)`, and `eta`, the shared inverse dispersion of that
#'     noise (only when `pathogen_noise = TRUE`).
#'   * `phi`: the inverse dispersion of the negative binomial observation
#'     model.
#'   * `tau` or `Sigma`: the smoothing standard deviation(s) or covariance
#'     matrix, depending on the smoothing structure.
#'   * `day_of_week_simplex`: the day-of-week effect (only when
#'     `dow_effect = TRUE`).
#' @export
#'
#' @srrstats {G1.4} uses `Roxygen2` documentation
#' @srrstats {G1.3, BS1.3} Arguments `n_chain`, `n_warmup`, `n_iter`, `thin`, and
#'   `adapt_delta` control the computational process and are described clearly
#'   in the documentation.
#' @srrstats {BS1.3a, BS2.8} To use parameter values from a previous fit as
#'   starting points for a new fit, they can be passed via the \code{init}
#'   parameter. See \code{?rstan::sampling} for details.
#' @srrstats {BS2.7, BS2.11} Starting values can be specified with the
#'   `init` argument to `rstan::sampling()` specified optionally here.
#' @srrstats {BS2.9, BS2.10} Stan uses random seeds that are different per chain
#' @srrstats {BS2.12} Argument `verbose`, which defaults to `TRUE`, controls the
#'   verbosity of output, showing all messages, warnings, errors, and progress
#'   indicators.
#' @srrstats {BS2.13} When `verbose = FALSE`, messages and progress indicators
#'   are suppressed (via `refresh = 0` and `show_messages = FALSE`) while
#'   warnings and errors are retained.
#' @srrstats {BS2.14} Argument `suppress_warnings`, which defaults to `FALSE`,
#'   enables suppression of warnings when set to `TRUE`.
#' @srrstats {BS2.15} Errors from Stan sampling are caught via `tryCatch` and
#'   re-raised with additional context. Stan failures (e.g., initialization
#'   failures) that don't throw errors are detected and converted to errors.
#'   The error object contains the original error message and the constructed
#'   model for inspection.
#' @srrstats {BS5.1} Returned list is of class `EpiStrainDynamics.fit`, and
#'   includes the fit object of class `stanfit` and the pre-specified
#'   constructed model object of class `EpiStrainDynamics.model`. The
#'   model object itself includes the input data.
#' @srrstats {BS5.0} Return values include mcmc components, including
#'   starting value(s) and seed(s).
#' @srrstats {BS6.4} `rstan` built-in summary methods can be used on the
#'   fitted model object
#' @srrstats {BS6.0} The fitting function, `fit_model()`, returns a list of
#'   class `EpiStrainDynamics.fit` as well as a class for the family of pathogen
#'   structure used (`ps`, `rw`, `ps_single`, or `rw_single`). The first item
#'   in the list is the stan output of class `stanfit`, which can be printed
#'   using default `rstan` print methods. By default the list is output to the
#'   console. A user can further interrogate the fit itself, or calculate
#'   epidemiological metrics using the `metrics` family of functions.
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
#' # Suppress progress and messages but keep warnings/errors
#' fit <- fit_model(mod, verbose = FALSE)
#'
#' # Suppress warnings too
#' fit <- fit_model(mod, verbose = FALSE, suppress_warnings = TRUE)
#'
#' # Catch errors and inspect
#' result <- tryCatch(
#'   fit_model(mod),
#'   error = function(e) e
#' )
#' if (inherits(result, "EpiStrainDynamics.fit.error")) {
#'   cat("Fitting failed:", result$message, "\n")
#'   # Can still access the model: result$constructed_model
#' }
fit_model <- function(constructed_model,
                      n_chain = 4,
                      n_iter = 2000,
                      n_warmup = floor(n_iter / 2),
                      thin = 1,
                      adapt_delta = 0.9,
                      multi_cores = TRUE,
                      verbose = TRUE,
                      suppress_warnings = FALSE,
                      seed = NULL, ...) {
  # validate inputs
  validate_class_inherits(constructed_model, "EpiStrainDynamics.model")
  validate_n_chain(n_chain)
  validate_n_iter(n_iter)
  validate_n_warmup(n_warmup)
  validate_thin(thin)
  validate_adapt_delta(adapt_delta)
  validate_seed(seed)
  validate_verbose(verbose)
  validate_suppress_warnings(suppress_warnings)
  validate_multi_cores(multi_cores)
  validate_mcmc_params_collective(n_iter, n_warmup, n_chain, thin, seed,
    suppress_warnings = suppress_warnings
  )

  UseMethod("fit_model")
}

#' Raise a catchable fitting error
#'
#' Signals an `EpiStrainDynamics.fit.error` condition carrying the constructed
#' model, so callers can catch a failed fit and inspect what was being fit.
#'
#' @noRd
abort_fit_error <- function(message, call, constructed_model, error_type) {
  cli::cli_abort(
    message,
    class = "EpiStrainDynamics.fit.error",
    call = call,
    constructed_model = constructed_model,
    error_type = error_type
  )
}

#' @rdname fit_model
#' @export
fit_model.rw_subtyped <- function(constructed_model,
                                  n_chain = 4,
                                  n_iter = 2000,
                                  n_warmup = floor(n_iter / 2),
                                  thin = 1,
                                  adapt_delta = 0.9,
                                  multi_cores = TRUE,
                                  verbose = TRUE,
                                  suppress_warnings = FALSE,
                                  seed = NULL, ...) {
  # Prepare Stan sampling call
  stan_call <- function() {
    rstan::sampling(
      stanmodels$rw_subtyped,
      data = constructed_model$standata,
      #' @srrstats {G2.4, G2.4a} explicit conversion to integers
      chains = as.integer(n_chain),
      iter = as.integer(n_iter),
      warmup = as.integer(n_warmup),
      thin = as.integer(thin),
      seed = ifelse(!is.null(seed), as.integer(seed),
        sample.int(.Machine$integer.max, 1)
      ),
      cores = ifelse(multi_cores == TRUE, parallel::detectCores() - 1, 1),
      control = list(adapt_delta = adapt_delta),
      #' @srrstats {BS2.12, BS2.13} refresh controls progress display
      refresh = ifelse(verbose == TRUE, 500, 0),
      #' @srrstats {BS2.13} show_messages controls informational messages
      show_messages = verbose,
      ...
    )
  }

  #' @srrstats {BS2.15} Catch errors from Stan and provide informative handling
  fit_object <- tryCatch(
    {
      #' @srrstats {BS2.14} Conditionally suppress warnings
      if (suppress_warnings) {
        suppressWarnings(stan_call())
      } else {
        stan_call()
      }
    },
    error = function(e) {
      # Create error object that can be inspected
      abort_fit_error(
        conditionMessage(e),
        conditionCall(e), constructed_model, class(e)
      )
    }
  )

  #' @srrstats {BS2.15} Check if Stan sampling actually succeeded
  #' Stan sometimes returns a failed stanfit object rather than throwing an error
  if (!inherits(fit_object, "stanfit") ||
    length(fit_object@sim$samples) == 0) {
    abort_fit_error(
      "Stan sampling failed to produce valid samples. Check for initialization failures or data issues.",
      sys.call(), constructed_model, "stan_failure"
    )
  }

  #' @srrstats {BS5.5} the `model` return object is of class `stanfit`, which
  #'   includes information about convergence
  out <- list(
    fit = fit_object,
    constructed_model = constructed_model,
    mcmc_info = list(
      seed = seed,
      chain_seeds = vapply(fit_object@stan_args, function(x) x$seed, integer(1)),
      init = rstan::get_inits(fit_object),
      n_chain = n_chain,
      n_iter = n_iter,
      n_warmup = n_warmup,
      thin = thin
    )
  )

  class(out) <- c("rw", "EpiStrainDynamics.fit", class(out))
  out
}

#' @rdname fit_model
#' @export
fit_model.ps_subtyped <- function(constructed_model,
                                  n_chain = 4,
                                  n_iter = 2000,
                                  n_warmup = floor(n_iter / 2),
                                  thin = 1,
                                  adapt_delta = 0.9,
                                  multi_cores = TRUE,
                                  verbose = TRUE,
                                  suppress_warnings = FALSE,
                                  seed = NULL, ...) {
  stan_call <- function() {
    rstan::sampling(
      stanmodels$ps_subtyped,
      data = constructed_model$standata,
      chains = as.integer(n_chain),
      iter = as.integer(n_iter),
      warmup = as.integer(n_warmup),
      thin = as.integer(thin),
      seed = ifelse(!is.null(seed), as.integer(seed),
        sample.int(.Machine$integer.max, 1)
      ),
      cores = ifelse(multi_cores == TRUE, parallel::detectCores() - 1, 1),
      control = list(adapt_delta = adapt_delta),
      refresh = ifelse(verbose == TRUE, 500, 0),
      show_messages = verbose,
      ...
    )
  }

  fit_object <- tryCatch(
    {
      if (suppress_warnings) {
        suppressWarnings(stan_call())
      } else {
        stan_call()
      }
    },
    error = function(e) {
      abort_fit_error(
        conditionMessage(e),
        conditionCall(e), constructed_model, class(e)
      )
    }
  )

  if (!inherits(fit_object, "stanfit") ||
    length(fit_object@sim$samples) == 0) {
    abort_fit_error(
      "Stan sampling failed to produce valid samples. Check for initialization failures or data issues.",
      sys.call(), constructed_model, "stan_failure"
    )
  }

  out <- list(
    fit = fit_object,
    constructed_model = constructed_model,
    mcmc_info = list(
      seed = seed,
      chain_seeds = vapply(fit_object@stan_args, function(x) x$seed, integer(1)),
      init = rstan::get_inits(fit_object),
      n_chain = n_chain,
      n_iter = n_iter,
      n_warmup = n_warmup,
      thin = thin
    )
  )

  class(out) <- c("ps", "EpiStrainDynamics.fit", class(out))

  out
}

#' @rdname fit_model
#' @export
fit_model.rw_multiple <- function(constructed_model,
                                  n_chain = 4,
                                  n_iter = 2000,
                                  n_warmup = floor(n_iter / 2),
                                  thin = 1,
                                  adapt_delta = 0.9,
                                  multi_cores = TRUE,
                                  verbose = TRUE,
                                  suppress_warnings = FALSE,
                                  seed = NULL, ...) {
  stan_call <- function() {
    rstan::sampling(
      stanmodels$rw_multiple,
      data = constructed_model$standata,
      chains = as.integer(n_chain),
      iter = as.integer(n_iter),
      warmup = as.integer(n_warmup),
      thin = as.integer(thin),
      seed = ifelse(!is.null(seed), as.integer(seed),
        sample.int(.Machine$integer.max, 1)
      ),
      cores = ifelse(multi_cores == TRUE, parallel::detectCores() - 1, 1),
      control = list(adapt_delta = adapt_delta),
      refresh = ifelse(verbose == TRUE, 500, 0),
      show_messages = verbose,
      ...
    )
  }

  fit_object <- tryCatch(
    {
      if (suppress_warnings) {
        suppressWarnings(stan_call())
      } else {
        stan_call()
      }
    },
    error = function(e) {
      abort_fit_error(
        conditionMessage(e),
        conditionCall(e), constructed_model, class(e)
      )
    }
  )

  if (!inherits(fit_object, "stanfit") ||
    length(fit_object@sim$samples) == 0) {
    abort_fit_error(
      "Stan sampling failed to produce valid samples. Check for initialization failures or data issues.",
      sys.call(), constructed_model, "stan_failure"
    )
  }

  out <- list(
    fit = fit_object,
    constructed_model = constructed_model,
    mcmc_info = list(
      seed = seed,
      chain_seeds = vapply(fit_object@stan_args, function(x) x$seed, integer(1)),
      init = rstan::get_inits(fit_object),
      n_chain = n_chain,
      n_iter = n_iter,
      n_warmup = n_warmup,
      thin = thin
    )
  )

  class(out) <- c("rw", "EpiStrainDynamics.fit", class(out))
  out
}

#' @rdname fit_model
#' @export
fit_model.ps_multiple <- function(constructed_model,
                                  n_chain = 4,
                                  n_iter = 2000,
                                  n_warmup = floor(n_iter / 2),
                                  thin = 1,
                                  adapt_delta = 0.9,
                                  multi_cores = TRUE,
                                  verbose = TRUE,
                                  suppress_warnings = FALSE,
                                  seed = NULL, ...) {
  stan_call <- function() {
    rstan::sampling(
      stanmodels$ps_multiple,
      data = constructed_model$standata,
      chains = as.integer(n_chain),
      iter = as.integer(n_iter),
      warmup = as.integer(n_warmup),
      thin = as.integer(thin),
      seed = ifelse(!is.null(seed), as.integer(seed),
        sample.int(.Machine$integer.max, 1)
      ),
      cores = ifelse(multi_cores == TRUE, parallel::detectCores() - 1, 1),
      control = list(adapt_delta = adapt_delta),
      refresh = ifelse(verbose == TRUE, 500, 0),
      show_messages = verbose,
      ...
    )
  }

  fit_object <- tryCatch(
    {
      if (suppress_warnings) {
        suppressWarnings(stan_call())
      } else {
        stan_call()
      }
    },
    error = function(e) {
      abort_fit_error(
        conditionMessage(e),
        conditionCall(e), constructed_model, class(e)
      )
    }
  )

  if (!inherits(fit_object, "stanfit") ||
    length(fit_object@sim$samples) == 0) {
    abort_fit_error(
      "Stan sampling failed to produce valid samples. Check for initialization failures or data issues.",
      sys.call(), constructed_model, "stan_failure"
    )
  }

  out <- list(
    fit = fit_object,
    constructed_model = constructed_model,
    mcmc_info = list(
      seed = seed,
      chain_seeds = vapply(fit_object@stan_args, function(x) x$seed, integer(1)),
      init = rstan::get_inits(fit_object),
      n_chain = n_chain,
      n_iter = n_iter,
      n_warmup = n_warmup,
      thin = thin
    )
  )

  class(out) <- c("ps", "EpiStrainDynamics.fit", class(out))
  out
}

#' @rdname fit_model
#' @export
fit_model.rw_single <- function(constructed_model,
                                n_chain = 4,
                                n_iter = 2000,
                                n_warmup = floor(n_iter / 2),
                                thin = 1,
                                adapt_delta = 0.9,
                                multi_cores = TRUE,
                                verbose = TRUE,
                                suppress_warnings = FALSE,
                                seed = NULL, ...) {
  stan_call <- function() {
    rstan::sampling(
      stanmodels$rw_single,
      data = constructed_model$standata,
      chains = as.integer(n_chain),
      iter = as.integer(n_iter),
      warmup = as.integer(n_warmup),
      thin = as.integer(thin),
      seed = ifelse(!is.null(seed), as.integer(seed),
        sample.int(.Machine$integer.max, 1)
      ),
      cores = ifelse(multi_cores == TRUE, parallel::detectCores() - 1, 1),
      control = list(adapt_delta = adapt_delta),
      refresh = ifelse(verbose == TRUE, 500, 0),
      show_messages = verbose,
      ...
    )
  }

  fit_object <- tryCatch(
    {
      if (suppress_warnings) {
        suppressWarnings(stan_call())
      } else {
        stan_call()
      }
    },
    error = function(e) {
      abort_fit_error(
        conditionMessage(e),
        conditionCall(e), constructed_model, class(e)
      )
    }
  )

  if (!inherits(fit_object, "stanfit") ||
    length(fit_object@sim$samples) == 0) {
    abort_fit_error(
      "Stan sampling failed to produce valid samples. Check for initialization failures or data issues.",
      sys.call(), constructed_model, "stan_failure"
    )
  }

  out <- list(
    fit = fit_object,
    constructed_model = constructed_model,
    mcmc_info = list(
      seed = seed,
      chain_seeds = vapply(fit_object@stan_args, function(x) x$seed, integer(1)),
      init = rstan::get_inits(fit_object),
      n_chain = n_chain,
      n_iter = n_iter,
      n_warmup = n_warmup,
      thin = thin
    )
  )

  class(out) <- c("rw_single", "EpiStrainDynamics.fit", class(out))
  out
}

#' @rdname fit_model
#' @export
fit_model.ps_single <- function(constructed_model,
                                n_chain = 4,
                                n_iter = 2000,
                                n_warmup = floor(n_iter / 2),
                                thin = 1,
                                adapt_delta = 0.9,
                                multi_cores = TRUE,
                                verbose = TRUE,
                                suppress_warnings = FALSE,
                                seed = NULL, ...) {
  stan_call <- function() {
    rstan::sampling(
      stanmodels$ps_single,
      data = constructed_model$standata,
      chains = as.integer(n_chain),
      iter = as.integer(n_iter),
      warmup = as.integer(n_warmup),
      thin = as.integer(thin),
      seed = ifelse(!is.null(seed), as.integer(seed),
        sample.int(.Machine$integer.max, 1)
      ),
      cores = ifelse(multi_cores == TRUE, parallel::detectCores() - 1, 1),
      control = list(adapt_delta = adapt_delta),
      refresh = ifelse(verbose == TRUE, 500, 0),
      show_messages = verbose,
      ...
    )
  }

  fit_object <- tryCatch(
    {
      if (suppress_warnings) {
        suppressWarnings(stan_call())
      } else {
        stan_call()
      }
    },
    error = function(e) {
      abort_fit_error(
        conditionMessage(e),
        conditionCall(e), constructed_model, class(e)
      )
    }
  )

  if (!inherits(fit_object, "stanfit") ||
    length(fit_object@sim$samples) == 0) {
    abort_fit_error(
      "Stan sampling failed to produce valid samples. Check for initialization failures or data issues.",
      sys.call(), constructed_model, "stan_failure"
    )
  }

  out <- list(
    fit = fit_object,
    constructed_model = constructed_model,
    mcmc_info = list(
      seed = seed,
      chain_seeds = vapply(fit_object@stan_args, function(x) x$seed, integer(1)),
      init = rstan::get_inits(fit_object),
      n_chain = n_chain,
      n_iter = n_iter,
      n_warmup = n_warmup,
      thin = thin
    )
  )

  class(out) <- c("ps_single", "EpiStrainDynamics.fit", class(out))
  out
}
