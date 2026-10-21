# Basic Units Support Tests

test_that("single() accepts and processes units objects ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5)
  )
  df$cases <- as_units(c(10, 20, 30, 40, 50), "count")

  # Should work without error
  result <- single(data = df, case_timeseries = "cases", time = "date")

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "single")

  # Units should be converted to numeric
  expect_type(result$data$case_timeseries, "double")
  expect_equal(as.numeric(result$data$case_timeseries), c(10, 20, 30, 40, 50))
})

test_that("multiple() handles units objects in case_timeseries ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    alpha = c(10, 20, 30, 40, 50),
    delta = c(5, 10, 15, 20, 25)
  )
  df$cases <- as_units(c(100, 200, 300, 400, 500), "count")

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_type(result$data$case_timeseries, "double")
  expect_equal(as.numeric(result$data$case_timeseries), c(100, 200, 300, 400, 500))
})

test_that("multiple() handles units objects in component pathogens ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(100, 200, 300, 400, 500)
  )
  df$alpha <- as_units(c(10, 20, 30, 40, 50), "count")
  df$delta <- as_units(c(5, 10, 15, 20, 25), "count")

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_type(result$data$component_pathogens, "double")

  # Check values were converted correctly
  expect_equal(result$data$component_pathogens[1, ], c(10, 20, 30, 40, 50))
  expect_equal(result$data$component_pathogens[2, ], c(5, 10, 15, 20, 25))
})

test_that("subtyped() handles units objects across all pathogen types ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    week = seq.Date(as.Date("2020-01-01"), by = "week", length.out = 5),
    inf_A = c(50, 100, 150, 200, 250)
  )
  df$ili <- as_units(c(100, 200, 300, 400, 500), "count")
  df$inf_H3N2 <- as_units(c(10, 20, 30, 40, 50), "count")
  df$inf_H1N1 <- as_units(c(15, 30, 45, 60, 75), "count")
  df$inf_B <- as_units(c(20, 40, 60, 80, 100), "count")
  df$other <- as_units(c(5, 10, 15, 20, 25), "count")

  result <- subtyped(
    data = df,
    case_timeseries = "ili",
    time = "week",
    parent_timeseries = "inf_A",
    subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
    other_pathogen_timeseries = c("inf_B", "other")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_type(result$data$case_timeseries, "double")
  expect_type(result$data$component_pathogens, "double")
  expect_type(result$data$subtyped, "double")
})

# Mixed Standard and Non-Standard Column Types
test_that("Functions handle mix of numeric and units columns ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(100, 200, 300, 400, 500), # regular numeric
    delta = c(5, 10, 15, 20, 25) # regular numeric
  )
  df$alpha <- as_units(c(10, 20, 30, 40, 50), "count") # units object

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  # Both should be converted to numeric
  expect_type(result$data$case_timeseries, "double")
  expect_type(result$data$component_pathogens, "double")
})

test_that("Functions handle all-numeric columns when units package loaded ", {
  skip_if_not_installed("units")
  library(units)

  # Even with units loaded, regular numeric columns should work fine
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(10, 20, 30, 40, 50)
  )

  result <- single(data = df, case_timeseries = "cases", time = "date")

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_type(result$data$case_timeseries, "double")
})

test_that("Units conversion preserves values correctly ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5)
  )

  # Create units with decimals
  df$cases <- as_units(c(10.5, 20.7, 30.2, 40.9, 50.1), "count")

  result <- single(data = df, case_timeseries = "cases", time = "date")

  expect_equal(
    as.numeric(result$data$case_timeseries),
    c(10.5, 20.7, 30.2, 40.9, 50.1),
    tolerance = 1e-10
  )
})

test_that("Units are properly stripped before modeling ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 50),
    cases = as_units(rpois(50, 50), "1"),
    alpha = as_units(rpois(50, 20), "1"),
    delta = as_units(rpois(50, 15), "1")
  )

  result <- multiple(
    data = df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  # The returned data should be plain numeric, not units objects
  expect_false(inherits(result$data$case_timeseries, "units"))
  expect_false(inherits(result$data$component_pathogens, "units"))
  expect_true(is.numeric(result$data$case_timeseries))
  expect_true(is.numeric(result$data$component_pathogens))
})

# Error Cases: Verify Appropriate Rejection
test_that("Functions appropriately reject unsupported column types ", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c("10", "20", "30", "40", "50") # character, not numeric
  )

  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "must be numeric or have units"
  )
})

test_that("Functions reject factor columns appropriately ", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = factor(c(10, 20, 30, 40, 50)) # factor, not numeric
  )

  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "must be numeric or have units"
  )
})

test_that("Functions reject logical columns appropriately ", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(TRUE, FALSE, TRUE, FALSE, TRUE) # logical, not numeric
  )

  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "must be numeric or have units"
  )
})

# Helper Function Tests

test_that("check_column_numeric accepts units objects ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    numeric_col = 1:10,
    units_col = as_units(1:10, "1"),
    char_col = as.character(1:10)
  )

  # Should pass for numeric
  expect_silent(check_column_numeric(df, "numeric_col"))

  # Should pass for units
  expect_silent(check_column_numeric(df, "units_col"))

  # Should error for character
  expect_error(
    check_column_numeric(df, "char_col"),
    "must be numeric or have units"
  )
})

# Integration with tsibble Validation

test_that("Units objects work with tsibble validation ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 10),
    cases = as_units(1:10, "1"),
    alpha = as_units(11:20, "1")
  )

  # Should pass all tsibble validations (no gaps, regular, ordered)
  result <- single(
    data = df,
    case_timeseries = "cases",
    time = "date"
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_length(result$data$case_timeseries, 10)
})

# Time Series Objects with Units

test_that("Functions handle time series objects with units columns ", {
  skip_if_not_installed("units")
  skip_if_not_installed("xts")
  library(units)
  library(xts)

  dates <- seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5)
  cases <- as_units(c(10, 20, 30, 40, 50), "count")

  ts_obj <- xts(data.frame(cases = cases), order.by = dates)

  result <- single(data = ts_obj, case_timeseries = "cases")

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_type(result$data$case_timeseries, "double")
})

# Edge Case: Units with NA

test_that("Functions catch NA values in units columns ", {
  skip_if_not_installed("units")
  library(units)

  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5)
  )
  df$cases <- as_units(c(10, 20, NA, 40, 50), "count")

  # Should error due to missing data, not due to units
  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "Missing values.*NA.*found"
  )
})

# Integer Type Handling

test_that("Functions handle integer columns (standard but non-double numeric) ", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = as.integer(c(10, 20, 30, 40, 50)) # integer, not double
  )

  result <- single(data = df, case_timeseries = "cases", time = "date")

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  # Should accept integers as they are numeric
  expect_true(is.numeric(result$data$case_timeseries))
})
