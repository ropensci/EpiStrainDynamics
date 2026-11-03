# Test suite for post-fit metric calculations
# Uses cached fitted models from setup-fitted-models.R

# ==============================================================================
# TESTS: Rt() FUNCTION
# ==============================================================================

test_that("Rt() validates inputs correctly", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # Wrong class
  expect_error(
    Rt(list(fit = "not a fit"), gi_dist = gi_simple, tau_max = 7),
    "Input must be of class"
  )

  # Invalid tau_max
  expect_error(
    Rt(fit_rw_single, gi_dist = gi_simple, tau_max = -1),
    "must be a positive number"
  )

  expect_error(
    Rt(fit_rw_single, gi_dist = gi_simple, tau_max = 1.5),
    "must be a whole number"
  )

  # Missing gi_dist
  expect_error(
    Rt(fit_rw_single, tau_max = 7),
    "missing"
  )

  # Invalid gi_dist
  expect_error(
    Rt(fit_rw_single, gi_dist = "not a function", tau_max = 7),
    "must be a function"
  )
})

test_that("Rt() works for single pathogen random walk model", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  result <- Rt(fit_rw_single, gi_dist = gi_simple, tau_max = 7)

  # Check structure
  expect_s3_class(result, "Rt")
  expect_s3_class(result, "EpiStrainDynamics.metric")
  expect_named(result, c("measure", "fit", "constructed_model"))

  # Check measure data frame
  expect_s3_class(result$measure, "data.frame")
  expect_named(
    result$measure,
    c(
      "time", "y", "lb_50", "ub_50", "lb_95",
      "ub_95", "prob_above", "pathogen"
    )
  )

  # Check dimensions - should start from time index = tau_max
  n_days <- nrow(fit_rw_single$constructed_model$validated_tsbl)
  expected_rows <- n_days - 7 + 1
  expect_equal(nrow(result$measure), expected_rows)

  # Check values are reasonable
  expect_true(all(result$measure$y > 0)) # Rt should be positive
  expect_true(all(result$measure$lb_95 > 0))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
  expect_true(all(result$measure$lb_50 >= result$measure$lb_95))

  # Check prob_above is between 0 and 1
  expect_true(all(result$measure$prob_above >= 0 & result$measure$prob_above <= 1))

  # Check pathogen name
  expect_equal(unique(result$measure$pathogen), "cases")
})

test_that("Rt() works for single pathogen p-spline model", {
  skip_if_not(exists("fit_ps_single"), "Cached fitted models not available")

  result <- Rt(fit_ps_single, gi_dist = gi_simple, tau_max = 7)

  expect_s3_class(result, "Rt")
  expect_s3_class(result$measure, "data.frame")

  n_days <- nrow(fit_ps_single$constructed_model$validated_tsbl)
  expected_rows <- n_days - 7 + 1
  expect_equal(nrow(result$measure), expected_rows)

  expect_true(all(result$measure$y > 0))
})

test_that("Rt() works for multiple pathogen models", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  result <- Rt(fit_rw_multi, gi_dist = gi_simple, tau_max = 7)

  expect_s3_class(result, "Rt")

  # Should have individual pathogens + Total
  pathogen_names <- unique(result$measure$pathogen)
  expect_true("Total" %in% pathogen_names)
  expect_true(all(c("alpha", "delta", "omicron", "other") %in% pathogen_names))

  n_days <- nrow(fit_rw_multi$constructed_model$validated_tsbl)
  expected_rows_per_pathogen <- n_days - 7 + 1
  expect_equal(nrow(result$measure), expected_rows_per_pathogen * 5)

  expect_true(all(result$measure$y > 0))
})

test_that("Rt() respects tau_max parameter", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  result_7 <- Rt(fit_rw_single, gi_dist = gi_simple, tau_max = 7)
  result_14 <- Rt(fit_rw_single, gi_dist = gi_simple, tau_max = 14)

  # Different tau_max should give different number of rows
  expect_lt(nrow(result_14$measure), nrow(result_7$measure))

  # Values should differ
  expect_false(all(abs(result_7$measure$y[1:10] - result_14$measure$y[1:10]) < 1e-10))
})

# ==============================================================================
# TESTS: incidence() FUNCTION
# ==============================================================================

test_that("incidence() validates inputs correctly", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  # Invalid class
  expect_error(
    incidence(list(fit = "not a fit"), dow = FALSE),
    "Input must be of class"
  )

  # dow = TRUE when model doesn't have dow_effect
  expect_error(
    incidence(fit_rw_single, dow = TRUE),
    "Day-of-week effects cannot be incorporated"
  )

  # dow must be NULL or logical
  expect_error(
    incidence(fit_rw_single, dow = "yes"),
    "must be.*NULL.*logical"
  )

  expect_error(
    incidence(fit_rw_single, dow = 1),
    "must be.*NULL.*logical"
  )
})

test_that("incidence() uses model's dow setting when dow = NULL", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  expect_message(
    result <- incidence(fit_rw_single, dow = NULL),
    "Using day-of-week setting from model.*FALSE"
  )
  expect_s3_class(result, "incidence")
})

test_that("incidence() default behavior uses model setting", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  expect_message(
    result <- incidence(fit_rw_single),
    "Using day-of-week setting from model"
  )
  expect_s3_class(result, "incidence")
})

test_that("incidence() explicit dow = FALSE works", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  expect_silent(
    result <- incidence(fit_rw_single, dow = FALSE)
  )
  expect_s3_class(result, "incidence")
})

test_that("incidence() works for single pathogen models", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  result <- incidence(fit_rw_single, dow = FALSE)

  expect_s3_class(result, "incidence")
  expect_s3_class(result, "EpiStrainDynamics.metric")

  expect_named(
    result$measure,
    c(
      "time", "y", "lb_50", "ub_50", "lb_95", "ub_95",
      "pathogen"
    )
  )

  n_days <- nrow(fit_rw_single$constructed_model$validated_tsbl)
  expect_equal(nrow(result$measure), n_days)

  expect_true(all(result$measure$y > 0))
  expect_true(all(result$measure$lb_95 > 0))

  # Check credible intervals are ordered
  expect_true(all(result$measure$lb_95 <= result$measure$lb_50))
  expect_true(all(result$measure$lb_50 <= result$measure$y))
  expect_true(all(result$measure$y <= result$measure$ub_50))
  expect_true(all(result$measure$ub_50 <= result$measure$ub_95))
})

test_that("incidence() works for multiple pathogen models", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  result <- incidence(fit_rw_multi, dow = FALSE)

  expect_s3_class(result, "incidence")

  pathogen_names <- unique(result$measure$pathogen)
  expect_true("Total" %in% pathogen_names)
  expect_equal(length(pathogen_names), 5)

  total_data <- result$measure[result$measure$pathogen == "Total", ]
  alpha_data <- result$measure[result$measure$pathogen == "alpha", ]
  expect_equal(nrow(total_data), nrow(alpha_data))
})

test_that("incidence() p-spline vs random walk gives different results", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")
  skip_if_not(exists("fit_ps_single"), "Cached fitted models not available")

  result_rw <- incidence(fit_rw_single, dow = FALSE)
  result_ps <- incidence(fit_ps_single, dow = FALSE)

  expect_equal(nrow(result_rw$measure), nrow(result_ps$measure))

  expect_false(all(abs(result_rw$measure$y - result_ps$measure$y) < 1e-10))
})

test_that("incidence() with dow = TRUE includes day-of-week effects", {
  skip_if_not(exists("fit_rw_single_dow"), "Model with dow_effect not available")

  result_with_dow <- incidence(fit_rw_single_dow, dow = TRUE)
  result_without_dow <- incidence(fit_rw_single_dow, dow = FALSE)

  expect_s3_class(result_with_dow, "incidence")
  expect_s3_class(result_without_dow, "incidence")

  # Results should differ when dow effects are included vs excluded
  expect_false(
    all(abs(result_with_dow$measure$y - result_without_dow$measure$y) < 1e-10)
  )
})

# ==============================================================================
# TESTS: growth_rate() FUNCTION
# ==============================================================================

test_that("growth_rate() validates inputs correctly", {
  expect_error(
    growth_rate(list(fit = "not a fit")),
    "Input must be of class"
  )
})

test_that("growth_rate() works for single pathogen models", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  result <- growth_rate(fit_rw_single)

  expect_s3_class(result, "growth_rate")
  expect_s3_class(result, "EpiStrainDynamics.metric")

  n_days <- nrow(fit_rw_single$constructed_model$validated_tsbl)
  expected_rows <- n_days - 2 + 1
  expect_equal(nrow(result$measure), expected_rows)

  # Growth rate can be positive or negative
  expect_true(any(result$measure$y > 0))
  expect_true(any(result$measure$y < 0))

  # But should be reasonable
  expect_true(all(abs(result$measure$y) < 10))
})

test_that("growth_rate() works for multiple pathogen models", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  result <- growth_rate(fit_rw_multi)

  expect_s3_class(result, "growth_rate")

  pathogen_names <- unique(result$measure$pathogen)
  expect_true("Total" %in% pathogen_names)
  expect_equal(length(pathogen_names), 5)
})

# ==============================================================================
# TESTS: proportion() FUNCTION
# ==============================================================================

test_that("proportion() validates inputs correctly", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  expect_error(
    proportion(list(fit = "not a fit")),
    "Input must be of class"
  )

  # Can't use proportion on single pathogen model
  expect_error(
    proportion(fit_rw_single),
    "Input must inherit from"
  )
})

test_that("proportion() works with default (all individual pathogens)", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  result <- proportion(fit_rw_multi)

  expect_s3_class(result, "proportion")
  expect_s3_class(result, "EpiStrainDynamics.metric")

  n_days <- nrow(fit_rw_multi$constructed_model$validated_tsbl)
  n_pathogens <- 4
  expected_rows <- n_days * n_pathogens
  expect_equal(nrow(result$measure), expected_rows)

  # Proportions should be between 0 and 1
  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))

  expect_true(all(result$measure$pathogen %in% c("alpha", "delta", "omicron", "other")))

  # For each time point, proportions should sum to ~1
  first_time <- result$measure[result$measure$time == result$measure$time[1], ]
  expect_true(abs(sum(first_time$y) - 1) < 0.1)
})

test_that("proportion() works with custom numerator and denominator", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  result <- proportion(fit_rw_multi,
    numerator_combination = "alpha",
    denominator_combination = c("alpha", "delta")
  )

  expect_s3_class(result, "proportion")

  n_days <- nrow(fit_rw_multi$constructed_model$validated_tsbl)
  expect_equal(nrow(result$measure), n_days)

  expect_true(all(result$measure$y >= 0))
  expect_true(all(result$measure$y <= 1))
})

test_that("proportion() records the resolved numerator/denominator combinations", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  default_result <- proportion(fit_rw_multi)
  expect_setequal(
    default_result$numerator_combination,
    c("alpha", "delta", "omicron", "other")
  )
  expect_setequal(
    default_result$denominator_combination,
    c("alpha", "delta", "omicron", "other")
  )

  custom_result <- proportion(fit_rw_multi,
    numerator_combination = "alpha",
    denominator_combination = c("alpha", "delta")
  )
  expect_equal(custom_result$numerator_combination, "alpha")
  expect_setequal(custom_result$denominator_combination, c("alpha", "delta"))
})

test_that("proportion() validates pathogen names", {
  skip_if_not(exists("fit_rw_multi"), "Cached fitted models not available")

  expect_error(
    proportion(fit_rw_multi,
      numerator_combination = "nonexistent",
      denominator_combination = "all"
    ),
    "invalid pathogen name"
  )

  expect_error(
    proportion(fit_rw_multi,
      numerator_combination = "alpha",
      denominator_combination = c("alpha", "nonexistent")
    ),
    "invalid pathogen name"
  )
})

# ==============================================================================
# TESTS: EDGE CASES AND ERROR HANDLING
# ==============================================================================

test_that("metrics handle missing values appropriately", {
  skip_if_not(exists("fit_rw_single"), "Cached fitted models not available")

  result <- incidence(fit_rw_single, dow = FALSE)

  # No NAs should be in output
  expect_false(any(is.na(result$measure$y)))
  expect_false(any(is.na(result$measure$lb_95)))
  expect_false(any(is.na(result$measure$ub_95)))
})

test_that("calc_stats() omits prob_above when threshold is NULL", {
  values <- c(-1, 0.5, 2, 3)

  expect_equal(calc_stats(values, threshold = 0)$prob_above, 0.75)
  expect_false("prob_above" %in% names(calc_stats(values, threshold = NULL)))
})

test_that("prob_above is returned for growth_rate() and Rt() only", {
  skip_if_not(exists("fit_rw_subtyped"), "Cached fitted models not available")

  # Incidence and proportions are always positive, so P(> 0) is always 1
  expect_false("prob_above" %in% names(incidence(fit_rw_subtyped)$measure))
  expect_false("prob_above" %in% names(proportion(fit_rw_subtyped)$measure))

  expect_true("prob_above" %in% names(growth_rate(fit_rw_subtyped)$measure))
  expect_true(
    "prob_above" %in% names(Rt(fit_rw_subtyped, gi_dist = gi_simple, tau_max = 7)$measure)
  )
})
