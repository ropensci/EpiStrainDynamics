# Test tsbox compatibility for EpiStrainDynamics intake functions
# Test that single(), multiple(), and subtyped() can accept any data.frame-like
# object handled by tsbox package, with optional time argument

# Helper function to create test data in different formats
create_test_ts_data <- function(class_type = "data.frame") {
  # Base data frame
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 100),
    cases = rpois(100, lambda = 50),
    alpha = rpois(100, lambda = 20),
    delta = rpois(100, lambda = 15),
    omicron = rpois(100, lambda = 10),
    other = rpois(100, lambda = 5)
  )

  # Convert to requested class type
  switch(class_type,
    "data.frame" = df,
    "data.table" = {
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::as.data.table(df)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "tibble" = {
      if (requireNamespace("tibble", quietly = TRUE)) {
        tibble::as_tibble(df)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "ts" = {
      # Convert to ts format (univariate time series)
      # For ts, we can only have one series, so return just cases
      ts(df$cases, start = 1, frequency = 1)
    },
    "mts" = {
      # Multivariate time series
      ts(df[, -1], start = 1, frequency = 1)
    },
    "xts" = {
      if (requireNamespace("xts", quietly = TRUE)) {
        xts::xts(df[, -1], order.by = df$date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "zoo" = {
      if (requireNamespace("zoo", quietly = TRUE)) {
        zoo::zoo(df[, -1], order.by = df$date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "zooreg" = {
      if (requireNamespace("zoo", quietly = TRUE)) {
        zoo::zooreg(df[, -1], start = as.Date("2020-01-01"), frequency = 1)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "tsibble" = {
      if (requireNamespace("tsibble", quietly = TRUE)) {
        tsibble::as_tsibble(df, index = date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "tibbletime" = {
      if (requireNamespace("tibbletime", quietly = TRUE)) {
        tibbletime::as_tbl_time(df, index = date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    stop("Unknown class type")
  )
}

# Test classes supported by tsbox and tsibble
# Separate data.frame-like classes from time series classes
dataframe_classes <- c("data.frame", "data.table", "tibble")
timeseries_classes <- c("mts", "xts", "zoo", "zooreg", "tsibble", "tibbletime")

test_that("single() accepts time series objects with optional time argument", {
  for (class_type in timeseries_classes) {
    test_data <- create_test_ts_data(class_type)

    # Test that single() can accept this data type without time argument
    expect_no_error(
      single(
        data = test_data,
        case_timeseries = "cases"
      )
    )

    # Should also work with explicit time = NULL
    expect_no_error(
      single(
        data = test_data,
        case_timeseries = "cases",
        time = NULL
      )
    )
  }
})

test_that("single() requires time argument for data.frame-like objects", {
  for (class_type in dataframe_classes) {
    test_data <- create_test_ts_data(class_type)

    # Should error when time is not provided for data.frame-like objects
    expect_error(
      single(
        data = test_data,
        case_timeseries = "cases"
      ),
      "When `time` is not specified, data must be a time series class object",
      info = paste("single() should require time for", class_type)
    )
  }
})

test_that("single() returns expected structure with all classes", {
  # Test with data frame (time required)
  test_df <- create_test_ts_data("data.frame")
  result_df <- single(
    data = test_df,
    case_timeseries = "cases",
    time = "date"
  )

  expect_s3_class(result_df, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result_df$pathogen_structure, "single")
  expect_equal(result_df$pathogen_names, "cases")
  expect_type(result_df$data, "list")
  expect_named(result_df$data, c("case_timeseries"))
  expect_length(result_df$data$case_timeseries, 100)

  # Test with time series object (time optional)
  if (requireNamespace("xts", quietly = TRUE)) {
    test_xts <- create_test_ts_data("xts")
    result_xts <- single(
      data = test_xts,
      case_timeseries = "cases"
    )

    expect_s3_class(result_xts, "EpiStrainDynamics.pathogen_structure")
    expect_equal(result_xts$pathogen_structure, "single")
    expect_equal(result_xts$pathogen_names, "cases")
    expect_length(result_xts$data$case_timeseries, 100)
  }
})

test_that("multiple() accepts data.frame-like objects with time argument", {
  for (class_type in dataframe_classes) {
    test_data <- create_test_ts_data(class_type)

    # Test that multiple() can accept this data type with time argument
    expect_no_error(
      multiple(
        data = test_data,
        case_timeseries = "cases",
        time = "date",
        component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
      )
    )
  }
})

test_that("multiple() accepts time series objs with optional time argument", {
  for (class_type in timeseries_classes) {
    test_data <- create_test_ts_data(class_type)

    # Test that multiple() can accept this data type without time argument
    expect_no_error(
      multiple(
        data = test_data,
        case_timeseries = "cases",
        component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
      )
    )
  }
})

test_that("multiple() requires time argument for data.frame-like objects", {
  for (class_type in dataframe_classes) {
    test_data <- create_test_ts_data(class_type)

    # Should error when time is not provided for data.frame-like objects
    expect_error(
      multiple(
        data = test_data,
        case_timeseries = "cases",
        component_pathogen_timeseries = c("alpha", "delta")
      ),
      "When `time` is not specified, data must be a time series class object"
    )
  }
})

test_that("multiple() returns expected structure with all classes", {
  # Test with data frame (time required)
  test_df <- create_test_ts_data("data.frame")
  result_df <- multiple(
    data = test_df,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
  )

  expect_s3_class(result_df, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result_df$pathogen_structure, "multiple")
  expect_equal(
    result_df$pathogen_names,
    c("alpha", "delta", "omicron", "other")
  )
  expect_type(result_df$data, "list")
  expect_named(result_df$data, c("case_timeseries", "component_pathogens"))
  expect_length(result_df$data$case_timeseries, 100)
  expect_equal(dim(result_df$data$component_pathogens), c(4, 100))

  # Test with time series object (time optional)
  if (requireNamespace("xts", quietly = TRUE)) {
    test_xts <- create_test_ts_data("xts")
    result_xts <- multiple(
      data = test_xts,
      case_timeseries = "cases",
      component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
    )

    expect_s3_class(result_xts, "EpiStrainDynamics.pathogen_structure")
    expect_equal(result_xts$pathogen_structure, "multiple")
    expect_equal(
      result_xts$pathogen_names,
      c("alpha", "delta", "omicron", "other")
    )
    expect_equal(dim(result_xts$data$component_pathogens), c(4, 100))
  }
})

# Helper to create subtyped test data
create_subtyped_test_data <- function(class_type = "data.frame") {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 100),
    cases = rpois(100, lambda = 50),
    inf_A = rpois(100, lambda = 20),
    inf_H1N1 = rpois(100, lambda = 10),
    inf_H3N2 = rpois(100, lambda = 8),
    inf_B = rpois(100, lambda = 15),
    other = rpois(100, lambda = 5)
  )

  # Convert to requested class type
  switch(class_type,
    "data.frame" = df,
    "data.table" = {
      if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::as.data.table(df)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "tibble" = {
      if (requireNamespace("tibble", quietly = TRUE)) {
        tibble::as_tibble(df)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "mts" = ts(df[, -1], start = 1, frequency = 1),
    "xts" = {
      if (requireNamespace("xts", quietly = TRUE)) {
        xts::xts(df[, -1], order.by = df$date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "zoo" = {
      if (requireNamespace("zoo", quietly = TRUE)) {
        zoo::zoo(df[, -1], order.by = df$date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "zooreg" = {
      if (requireNamespace("zoo", quietly = TRUE)) {
        zoo::zooreg(df[, -1], start = as.Date("2020-01-01"), frequency = 1)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "tsibble" = {
      if (requireNamespace("tsibble", quietly = TRUE)) {
        tsibble::as_tsibble(df, index = date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    "tibbletime" = {
      if (requireNamespace("tibbletime", quietly = TRUE)) {
        tibbletime::as_tbl_time(df, index = date)
      } else {
        skip(paste(class_type, "package not available"))
      }
    },
    stop("Unknown class type")
  )
}

test_that("subtyped() accepts data.frame-like objects with time argument", {
  for (class_type in dataframe_classes) {
    test_data <- create_subtyped_test_data(class_type)

    # Test that subtyped() can accept this data type with time argument
    expect_no_error(
      subtyped(
        data = test_data,
        case_timeseries = "cases",
        time = "date",
        parent_timeseries = "inf_A",
        subtyped_timeseries = c("inf_H1N1", "inf_H3N2"),
        other_pathogen_timeseries = c("inf_B", "other")
      )
    )
  }
})

test_that("subtyped() accepts time series objects with optional time arg", {
  for (class_type in timeseries_classes) {
    test_data <- create_subtyped_test_data(class_type)

    # Test that subtyped() can accept this data type without time argument
    expect_no_error(
      subtyped(
        data = test_data,
        case_timeseries = "cases",
        parent_timeseries = "inf_A",
        subtyped_timeseries = c("inf_H1N1", "inf_H3N2"),
        other_pathogen_timeseries = c("inf_B", "other")
      )
    )
  }
})

test_that("subtyped() requires time argument for data.frame-like objects", {
  for (class_type in dataframe_classes) {
    test_data <- create_subtyped_test_data(class_type)

    # Should error when time is not provided for data.frame-like objects
    expect_error(
      subtyped(
        data = test_data,
        case_timeseries = "cases",
        parent_timeseries = "inf_A",
        subtyped_timeseries = c("inf_H1N1", "inf_H3N2"),
        other_pathogen_timeseries = c("inf_B", "other")
      ),
      "When `time` is not specified, data must be a time series class object"
    )
  }
})

test_that("subtyped() returns expected structure with all classes", {
  # Test with data frame (time required)
  test_df <- create_subtyped_test_data("data.frame")
  result_df <- subtyped(
    data = test_df,
    case_timeseries = "cases",
    time = "date",
    parent_timeseries = "inf_A",
    subtyped_timeseries = c("inf_H1N1", "inf_H3N2"),
    other_pathogen_timeseries = c("inf_B", "other")
  )

  expect_s3_class(result_df, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result_df$pathogen_structure, "subtyped")
  expect_equal(
    result_df$pathogen_names,
    c("inf_H1N1", "inf_H3N2", "inf_B", "other")
  )
  expect_type(result_df$data, "list")
  expect_named(result_df$data, c("case_timeseries", "component_pathogens", "subtyped"))
  expect_length(result_df$data$case_timeseries, 100)
  expect_equal(dim(result_df$data$component_pathogens), c(3, 100))
  expect_equal(dim(result_df$data$subtyped), c(2, 100))

  # Test with time series object (time optional)
  if (requireNamespace("xts", quietly = TRUE)) {
    test_xts <- create_subtyped_test_data("xts")
    result_xts <- subtyped(
      data = test_xts,
      case_timeseries = "cases",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H1N1", "inf_H3N2"),
      other_pathogen_timeseries = c("inf_B", "other")
    )

    expect_s3_class(result_xts, "EpiStrainDynamics.pathogen_structure")
    expect_equal(result_xts$pathogen_structure, "subtyped")
    expect_equal(
      result_xts$pathogen_names,
      c("inf_H1N1", "inf_H3N2", "inf_B", "other")
    )
    expect_equal(dim(result_xts$data$component_pathogens), c(3, 100))
    expect_equal(dim(result_xts$data$subtyped), c(2, 100))
  }
})

# Integration tests with construct_model()

test_that("intake functions work with construct_model()", {
  # Test a subset of classes for integration (to keep tests fast)
  integration_classes <- c("data.frame", "tibble")

  for (class_type in integration_classes) {
    if (class_type == "tibble" && !requireNamespace("tibble", quietly = TRUE)) {
      skip("tibble not available")
    }

    test_data <- create_test_ts_data(class_type)

    # Test that construct_model works with the intake function output
    expect_no_error(
      construct_model(
        pathogen_structure = multiple(
          data = test_data,
          case_timeseries = "cases",
          time = "date",
          component_pathogen_timeseries = c(
            "alpha", "delta",
            "omicron", "other"
          )
        ),
        method = random_walk(),
        smoothing_params = smoothing_structure(
          "independent",
          tau_mean = rep(0, 4),
          tau_sd = rep(1, 4)
        ),
        dispersion_params = dispersion_structure(phi_mean = 0, phi_sd = 1),
        pathogen_noise = FALSE,
        dow_effect = TRUE
      )
    )
  }

  # Test with time series object (time optional)
  if (requireNamespace("xts", quietly = TRUE)) {
    test_xts <- create_test_ts_data("xts")

    expect_no_error(
      construct_model(
        pathogen_structure = multiple(
          data = test_xts,
          case_timeseries = "cases",
          component_pathogen_timeseries = c(
            "alpha", "delta",
            "omicron", "other"
          )
        ),
        method = random_walk(),
        smoothing_params = smoothing_structure(
          "independent",
          tau_mean = rep(0, 4),
          tau_sd = rep(1, 4)
        ),
        dispersion_params = dispersion_structure(phi_mean = 0, phi_sd = 1),
        pathogen_noise = FALSE,
        dow_effect = TRUE
      )
    )
  }
})

# Data integrity tests

test_that("intake functions preserve data integrity across conversions", {
  # Create original data
  original_df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 50),
    cases = rpois(50, lambda = 50),
    alpha = rpois(50, lambda = 20)
  )

  # Convert to different formats and check values are preserved
  for (class_type in c("data.frame", "tibble", "data.table")) {
    if (class_type %in% c("tibble", "data.table")) {
      pkg <- ifelse(class_type == "tibble", "tibble", "data.table")
      if (!requireNamespace(pkg, quietly = TRUE)) {
        skip(paste(pkg, "not available"))
      }
    }

    # Convert data
    if (class_type == "tibble") {
      test_data <- tibble::as_tibble(original_df)
    } else if (class_type == "data.table") {
      test_data <- data.table::as.data.table(original_df)
    } else {
      test_data <- original_df
    }

    # Use intake function
    result <- single(
      data = test_data,
      case_timeseries = "cases",
      time = "date"
    )

    # Check values are preserved
    expect_equal(
      result$data$case_timeseries,
      original_df$cases,
      info = paste("Values should be preserved for", class_type)
    )
  }

  # Test with time series object
  if (requireNamespace("xts", quietly = TRUE)) {
    test_xts <- xts::xts(original_df[, c("cases", "alpha")],
      order.by = original_df$date
    )

    result_xts <- single(
      data = test_xts,
      case_timeseries = "cases"
    )

    # Check values are preserved
    expect_equal(
      result_xts$data$case_timeseries,
      original_df$cases,
      info = "Values should be preserved for xts"
    )
  }
})

test_that("time series objects work correctly without time argument", {
  # Create test data
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 30),
    cases = rpois(30, lambda = 50),
    alpha = rpois(30, lambda = 20),
    delta = rpois(30, lambda = 15)
  )

  if (requireNamespace("xts", quietly = TRUE)) {
    # Convert to xts
    test_xts <- xts::xts(df[, -1], order.by = df$date)

    # Test single()
    result_single <- single(data = test_xts, case_timeseries = "cases")
    expect_length(result_single$data$case_timeseries, 30)

    # Test multiple()
    result_multiple <- multiple(
      data = test_xts,
      case_timeseries = "cases",
      component_pathogen_timeseries = c("alpha", "delta")
    )
    expect_equal(dim(result_multiple$data$component_pathogens), c(2, 30))
  }
})
