# Test helper function: check_list_columns
test_that("check_list_columns detects list columns in data.frame", {
  # Create data with list column
  df <- data.frame(
    time = 1:5,
    cases = c(10, 20, 30, 40, 50),
    metadata = I(list(
      list(a = 1), list(a = 2), list(a = 3), list(a = 4), list(a = 5)
    ))
  )

  expect_error(
    check_list_columns(df, c("time", "cases", "metadata")),
    "List columns detected"
  )

  expect_error(
    check_list_columns(df, c("time", "cases", "metadata")),
    "List columns are not supported for time series analysis"
  )
})

test_that("check_list_columns passes when no list columns present", {
  df <- data.frame(
    time = 1:5,
    cases = c(10, 20, 30, 40, 50),
    name = c("a", "b", "c", "d", "e")
  )

  expect_silent(
    check_list_columns(df, c("time", "cases", "name"))
  )
})

test_that("check_list_columns ignores non-relevant list columns", {
  df <- data.frame(
    time = 1:5,
    cases = c(10, 20, 30, 40, 50),
    metadata = I(list(
      list(a = 1), list(a = 2), list(a = 3), list(a = 4), list(a = 5)
    ))
  )

  # Should pass because metadata is not in relevant_cols
  expect_silent(
    check_list_columns(df, c("time", "cases"))
  )
})

test_that("check_list_columns handles multiple list columns", {
  df <- data.frame(
    time = 1:5,
    cases = c(10, 20, 30, 40, 50),
    metadata1 = I(list(
      list(a = 1), list(a = 2), list(a = 3), list(a = 4), list(a = 5)
    )),
    metadata2 = I(list(
      list(b = 1), list(b = 2), list(b = 3), list(b = 4), list(b = 5)
    ))
  )

  expect_error(
    check_list_columns(df, c("time", "cases", "metadata1", "metadata2")),
    "List columns detected"
  )
})

test_that("check_list_columns handles tibbles with list columns", {
  skip_if_not_installed("tibble")
  library(tibble)

  df <- tibble(
    time = 1:5,
    cases = c(10, 20, 30, 40, 50),
    metadata = list(
      list(a = 1), list(a = 2), list(a = 3), list(a = 4), list(a = 5)
    )
  )

  expect_error(
    check_list_columns(df, c("time", "cases", "metadata")),
    "List columns detected"
  )
})

# Test single() function with list columns
test_that("single() rejects data.frame with list columns", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(10, 20, 30, 40, 50),
    metadata = I(list(
      list(a = 1), list(a = 2), list(a = 3), list(a = 4), list(a = 5)
    ))
  )

  expect_error(
    single(data = df, case_timeseries = "metadata", time = "date"),
    "List columns detected"
  )
})

# Test multiple() function with list columns
test_that("multiple() rejects data.frame with list columns in
          component_pathogen_timeseries", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(100, 200, 300, 400, 500),
    alpha = c(10, 20, 30, 40, 50),
    delta = I(list(
      list(val = 5), list(val = 10), list(val = 15), list(val = 20), list(val = 25)
    ))
  )

  expect_error(
    multiple(
      data = df,
      case_timeseries = "cases",
      time = "date",
      component_pathogen_timeseries = c("alpha", "delta")
    ),
    "List columns detected"
  )
})

# Test subtyped() function with list columns
test_that("subtyped() rejects data.frame with list columns", {
  df <- data.frame(
    week = seq.Date(as.Date("2020-01-01"), by = "week", length.out = 5),
    ili = c(100, 200, 300, 400, 500),
    inf_A = c(50, 100, 150, 200, 250),
    inf_H3N2 = I(list(
      list(val = 10), list(val = 20), list(val = 30), list(val = 40), list(val = 50)
    )),
    inf_H1N1 = c(15, 30, 45, 60, 75),
    inf_B = c(20, 40, 60, 80, 100),
    other = c(5, 10, 15, 20, 25)
  )

  expect_error(
    subtyped(
      data = df,
      case_timeseries = "ili",
      time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
      other_pathogen_timeseries = c("inf_B", "other")
    ),
    "List columns detected"
  )
})

# Edge case: list column in time variable
test_that("Functions reject list column in time variable", {
  df <- data.frame(
    date = I(list(
      as.Date("2020-01-01"), as.Date("2020-01-02"),
      as.Date("2020-01-03"), as.Date("2020-01-04"), as.Date("2020-01-05")
    )),
    cases = c(10, 20, 30, 40, 50)
  )

  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "List columns detected"
  )
})

test_that("Functions handle list columns in data.table", {
  skip_if_not_installed("data.table")
  library(data.table)

  dt <- data.table(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = c(10, 20, 30, 40, 50),
    metadata = list(
      list(a = 1), list(a = 2), list(a = 3), list(a = 4), list(a = 5)
    )
  )

  expect_error(
    single(data = dt, case_timeseries = "metadata", time = "date"),
    "List columns detected"
  )
})

# Test informative error messages
test_that("Error message provides helpful guidance", {
  df <- data.frame(
    date = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 5),
    cases = I(list(
      c(10), c(20), c(30), c(40), c(50)
    ))
  )

  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "convert these columns to appropriate vector types"
  )

  expect_error(
    single(data = df, case_timeseries = "cases", time = "date"),
    "unlist()"
  )
})
