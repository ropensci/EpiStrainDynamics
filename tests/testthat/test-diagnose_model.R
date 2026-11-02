# Tests for diagnose_model() function
# Uses same fixtures as test-metrics.R to avoid duplication

# ==============================================================================
# TESTS: INPUT VALIDATION
# ==============================================================================

test_that("diagnose_model() validates input class", {
  expect_error(
    diagnose_model(list(fit = "not a fit")),
    "must be an .*EpiStrainDynamics.fit.* object"
  )

  expect_error(
    diagnose_model("not a model"),
    "must be an .*EpiStrainDynamics.fit.* object"
  )

  expect_error(
    diagnose_model(data.frame(x = 1, y = 2)),
    "must be an .*EpiStrainDynamics.fit.* object"
  )
})

test_that("diagnose_model() validates threshold parameters", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # rhat_threshold should be numeric
  expect_error(
    diagnose_model(fit_rw_single, rhat_threshold = "bad"),
    "must be numeric"
  )

  # eff_sample_threshold should be numeric
  expect_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = "bad"),
    "must be numeric"
  )

  # rhat_threshold must be >= 1.0
  expect_error(
    diagnose_model(fit_rw_single, rhat_threshold = 0.5),
    "must be greater than or equal to 1.0"
  )

  expect_error(
    diagnose_model(fit_rw_single, rhat_threshold = 0.99),
    "must be greater than or equal to 1.0"
  )

  # eff_sample_threshold must be positive
  expect_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = 0),
    "must be positive"
  )

  expect_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = -100),
    "must be positive"
  )

  # Should accept valid values
  expect_no_error(
    diagnose_model(fit_rw_single, rhat_threshold = 1.0)
  )

  expect_no_error(
    diagnose_model(fit_rw_single, rhat_threshold = 1.1)
  )

  expect_no_error(
    diagnose_model(fit_rw_single, eff_sample_threshold = 100)
  )
})

test_that("diagnose_model() warns about unusual threshold values", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # Very high rhat_threshold should warn
  expect_warning(
    diagnose_model(fit_rw_single, rhat_threshold = 5.0),
    "unusually high"
  )

  # Very low eff_sample_threshold should warn
  expect_warning(
    diagnose_model(fit_rw_single, eff_sample_threshold = 10),
    "very low"
  )

  # Very high eff_sample_threshold should warn
  expect_warning(
    diagnose_model(fit_rw_single, eff_sample_threshold = 10000),
    "very high"
  )
})

# ==============================================================================
# TESTS: OUTPUT STRUCTURE
# ==============================================================================

test_that("diagnose_model() returns correct structure", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  # Should return a list
  expect_type(diag, "list")

  # Should have expected components
  expect_named(diag, c(
    "convergence", "rhat_issues", "eff_sample_issues",
    "max_rhat", "min_neff", "summary", "parameters",
    "rhat_threshold", "eff_sample_threshold"
  ))
  expect_s3_class(diag, "EpiStrainDynamics.diagnostics")
})

test_that("diagnose_model() convergence is logical", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  expect_type(diag$convergence, "logical")
  expect_length(diag$convergence, 1)
  expect_false(is.na(diag$convergence))
})

test_that("diagnose_model() issues are character vectors", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  expect_type(diag$rhat_issues, "character")
  expect_type(diag$eff_sample_issues, "character")

  # Can be empty (no issues)
  expect_true(length(diag$rhat_issues) >= 0)
  expect_true(length(diag$eff_sample_issues) >= 0)
})

test_that("diagnose_model() max_rhat and min_neff are numeric", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  expect_type(diag$max_rhat, "double")
  expect_type(diag$min_neff, "double")

  expect_length(diag$max_rhat, 1)
  expect_length(diag$min_neff, 1)

  expect_true(is.finite(diag$max_rhat))
  expect_true(is.finite(diag$min_neff))
})

test_that("diagnose_model() summary is a matrix", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  expect_true(is.matrix(diag$summary))

  # Should have standard Stan summary columns
  expected_cols <- c(
    "mean", "se_mean", "sd", "2.5%", "25%", "50%",
    "75%", "97.5%", "n_eff", "Rhat"
  )
  expect_true(all(expected_cols %in% colnames(diag$summary)))
})

# ==============================================================================
# TESTS: CONVERGENCE LOGIC
# ==============================================================================

test_that("diagnose_model() convergence logic is correct", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single,
    rhat_threshold = 1.1,
    eff_sample_threshold = 100
  )

  # If convergence is TRUE, there should be no issues
  if (diag$convergence) {
    expect_equal(length(diag$rhat_issues), 0)
    expect_equal(length(diag$eff_sample_issues), 0)
  }

  # If convergence is FALSE, there should be at least one issue
  if (!diag$convergence) {
    total_issues <- length(diag$rhat_issues) + length(diag$eff_sample_issues)
    expect_gt(total_issues, 0)
  }
})

# ==============================================================================
# TESTS: WORKS WITH DIFFERENT MODEL TYPES
# ==============================================================================

test_that("diagnose_model() works with single pathogen random walk", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  expect_no_error(diag <- diagnose_model(fit_rw_single))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

test_that("diagnose_model() works with single pathogen p-spline", {
  skip_if_not(exists("fit_ps_single"), "Cached fitted models not available")

  expect_no_error(diag <- diagnose_model(fit_ps_single))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

test_that("diagnose_model() works with multiple pathogen random walk", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  expect_no_error(diag <- diagnose_model(fit_rw_multi))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

test_that("diagnose_model() works with multiple pathogen p-spline", {
  skip_if_not(exists("fit_ps_multi"), "Cached fitted models not available")

  expect_no_error(diag <- diagnose_model(fit_ps_multi))
  expect_type(diag, "list")
  expect_type(diag$convergence, "logical")
})

# ==============================================================================
# TESTS: DIAGNOSTIC VALUES ARE REASONABLE
# ==============================================================================

test_that("diagnose_model() max_rhat is >= 1", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  # R-hat should always be >= 1 by definition
  expect_gte(diag$max_rhat, 1.0)
})

test_that("diagnose_model() min_neff is positive", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  # Effective sample size should be positive
  expect_gt(diag$min_neff, 0)
})

test_that("diagnose_model() identifies parameter names in issues", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single,
    rhat_threshold = 1.01,
    eff_sample_threshold = 500
  )

  # If there are issues, they should be named parameters
  if (length(diag$rhat_issues) > 0) {
    expect_true(all(nzchar(diag$rhat_issues)))

    # Should be valid row names from summary
    expect_true(all(diag$rhat_issues %in% rownames(diag$summary)))
  }

  if (length(diag$eff_sample_issues) > 0) {
    expect_true(all(nzchar(diag$eff_sample_issues)))
    expect_true(all(diag$eff_sample_issues %in% rownames(diag$summary)))
  }
})

# ==============================================================================
# TESTS: CONSOLE OUTPUT (captured)
# ==============================================================================

test_that("diagnose_model() prints summary to console", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  output <- capture.output(diag <- diagnose_model(fit_rw_single))

  # Should have printed something
  expect_gt(length(output), 0)

  # Should contain key phrases
  expect_true(any(grepl("Model Convergence Diagnostics", output)))
  expect_true(any(grepl("Overall convergence", output)))
  expect_true(any(grepl("Maximum R-hat", output)))
  expect_true(any(grepl("Minimum n_eff", output)))
})

test_that("diagnose_model() shows GOOD or ISSUES DETECTED", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  output <- capture.output(diag <- diagnose_model(fit_rw_single))

  # Should contain either GOOD or ISSUES DETECTED
  output_text <- paste(output, collapse = " ")
  expect_true(grepl("GOOD|ISSUES DETECTED", output_text))
})

test_that("diagnose_model() console output matches convergence status", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  output <- capture.output(diag <- diagnose_model(fit_rw_single))
  output_text <- paste(output, collapse = " ")

  if (diag$convergence) {
    expect_true(grepl("GOOD", output_text))
  } else {
    expect_true(grepl("ISSUES DETECTED", output_text))
  }
})

# ==============================================================================
# TESTS: INVISIBLE RETURN
# ==============================================================================

test_that("diagnose_model() returns invisibly", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # The function should return invisible
  # We can't directly test invisibility, but we can check that
  # assigning doesn't print
  output <- capture.output({
    result <- diagnose_model(fit_rw_single)
  })

  # Should have printed diagnostics but not returned value
  expect_true(any(grepl("Model Convergence", output)))
})

# ==============================================================================
# TESTS: COMPARISON ACROSS MODEL TYPES
# ==============================================================================

test_that("diagnose_model() produces consistent output across model types", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag_rw_single <- diagnose_model(fit_rw_single)
  diag_ps_single <- diagnose_model(fit_ps_single)
  diag_rw_multi <- diagnose_model(fit_rw_multi)
  diag_ps_multi <- diagnose_model(fit_ps_multi)

  # All should have the same structure
  expected_names <- c(
    "convergence", "rhat_issues", "eff_sample_issues",
    "max_rhat", "min_neff", "summary", "parameters",
    "rhat_threshold", "eff_sample_threshold"
  )

  expect_named(diag_rw_single, expected_names)
  expect_named(diag_ps_single, expected_names)
  expect_named(diag_rw_multi, expected_names)
  expect_named(diag_ps_multi, expected_names)
})

test_that("diagnose_model() summary has more parameters for multi models", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag_single <- diagnose_model(fit_rw_single)
  diag_multi <- diagnose_model(fit_rw_multi)

  # Multi-pathogen model should have more parameters
  n_params_single <- nrow(diag_single$summary)
  n_params_multi <- nrow(diag_multi$summary)

  expect_gt(n_params_multi, n_params_single)
})

# ==============================================================================
# TESTS: COMPARISON WITH RSTAN DIAGNOSTICS
# ==============================================================================

test_that("diagnose_model() max_rhat matches Stan summary", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  # Get R-hat values directly from Stan
  stan_summary <- rstan::summary(fit_rw_single$fit)$summary
  stan_rhat <- stan_summary[, "Rhat"]
  stan_max_rhat <- max(stan_rhat, na.rm = TRUE)

  # Should match
  expect_equal(diag$max_rhat, stan_max_rhat)
})

test_that("diagnose_model() min_neff matches Stan summary", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  diag <- diagnose_model(fit_rw_single)

  # Get n_eff values directly from Stan
  stan_summary <- rstan::summary(fit_rw_single$fit)$summary
  stan_neff <- stan_summary[, "n_eff"]
  stan_min_neff <- min(stan_neff, na.rm = TRUE)

  # Should match
  expect_equal(diag$min_neff, stan_min_neff)
})

# ==============================================================================
# TESTS: CONVERTING POSITIONS TO DATA TIME
# ==============================================================================

test_that("index_to_time() steps Date and POSIXct times by the observation spacing", {
  weekly <- as.Date("2024-01-01") + 7 * 0:4
  expect_equal(
    index_to_time(c(-1, 1, 3, 6), weekly),
    as.Date("2024-01-01") + 7 * c(-2, 0, 2, 5)
  )

  hourly <- as.POSIXct("2024-01-01 00:00", tz = "UTC") + 3600 * 0:4
  expect_equal(
    index_to_time(c(1, 2.5), hourly),
    as.POSIXct("2024-01-01 00:00", tz = "UTC") + 3600 * c(0, 1.5)
  )
})

test_that("index_to_time() rounds yearweek, yearmonth and yearquarter times to whole units", {
  weeks <- tsibble::yearweek("2024 W01") + 0:4
  expect_no_warning(res <- index_to_time(c(-1, 1, 2.4, 6), weeks))
  expect_s3_class(res, "yearweek")
  expect_equal(res, tsibble::yearweek("2024 W01") + c(-2, 0, 1, 5))

  fortnights <- tsibble::yearweek("2024 W01") + 2 * 0:4
  expect_equal(
    index_to_time(c(1, 2.4), fortnights),
    tsibble::yearweek("2024 W01") + c(0, 3)
  )

  months <- tsibble::yearmonth("2024 Jan") + 0:4
  expect_no_warning(res <- index_to_time(c(0.4, 2.4), months))
  expect_equal(res, tsibble::yearmonth("2024 Jan") + c(-1, 1))

  quarters <- tsibble::yearquarter("2024 Q1") + 0:4
  expect_no_warning(res <- index_to_time(c(0.4, 2.4), quarters))
  expect_equal(res, tsibble::yearquarter("2024 Q1") + c(-1, 1))
})

test_that("index_to_time() keeps integer times as integers", {
  res <- index_to_time(c(-1.4, 2.4), 1:5)
  expect_type(res, "integer")
  expect_equal(res, c(-1L, 2L))

  expect_equal(index_to_time(c(1, 2.5), c(1, 2, 3)), c(1, 2.5))
})

test_that("index_to_time() returns NA with one warning when times cannot be converted", {
  expect_warning(
    res <- index_to_time(c(1, 2), c("a", "b")),
    "Could not convert spline coefficient positions"
  )
  expect_equal(res, c(NA_character_, NA_character_))

  expect_equal(index_to_time(c(1, 2), as.Date("2024-01-01")), as.Date(c(NA, NA)))
})
