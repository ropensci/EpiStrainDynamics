# Pathogen structure tests

# Tests for single() function
test_that("single() creates correct structure", {
  check_package_data()

  result <- single(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date"
  )

  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "single")
  expect_equal(result$pathogen_names, "cases")
  expect_equal(length(result$data$case_timeseries), nrow(sarscov2))
  expect_type(result$data$case_timeseries, "double")
})

test_that("single() validates column existence and types", {
  check_package_data()

  # Test non-existent time column
  expect_error(
    single(data = sarscov2, case_timeseries = "cases", time = "nonexistent"),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent case column
  expect_error(
    single(data = sarscov2, case_timeseries = "nonexistent", time = "date"),
    "Column .nonexistent. not found in data"
  )

  # Test non-numeric case column (if you have a non-numeric column)
  test_data <- sarscov2
  test_data$text_col <- as.character(test_data$cases)
  expect_error(
    single(data = test_data, case_timeseries = "text_col", time = "date"),
    "numeric"
  )
})

# Tests for multiple() function
test_that("multiple() creates correct structure", {
  check_package_data()

  result <- multiple(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta", "omicron", "other")
  )

  # Test basic structure
  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "multiple")
  expect_equal(result$pathogen_names, c("alpha", "delta", "omicron", "other"))

  # Test data structure
  expect_equal(length(result$data$case_timeseries), nrow(sarscov2))
  expect_equal(nrow(result$data$component_pathogens), 4)
  expect_equal(ncol(result$data$component_pathogens), nrow(sarscov2))

  # Verify matrix transposition is correct
  expect_equal(result$data$component_pathogens[1, ], sarscov2$alpha)
  expect_equal(result$data$component_pathogens[2, ], sarscov2$delta)
  expect_equal(result$data$component_pathogens[3, ], sarscov2$omicron)
  expect_equal(result$data$component_pathogens[4, ], sarscov2$other)
})

test_that("multiple() handles single component pathogen", {
  check_package_data()

  result <- multiple(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha")
  )

  expect_equal(result$pathogen_names, "alpha")
  expect_equal(nrow(result$data$component_pathogens), 1)
  expect_equal(ncol(result$data$component_pathogens), nrow(sarscov2))
  expect_equal(result$data$component_pathogens[1, ], sarscov2$alpha)
})

test_that("multiple() validates column existence and types", {
  check_package_data()

  # Test non-existent time column
  expect_error(
    multiple(
      data = sarscov2,
      case_timeseries = "cases",
      time = "nonexistent",
      component_pathogen_timeseries = c("alpha")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent case column
  expect_error(
    multiple(
      data = sarscov2,
      case_timeseries = "nonexistent",
      time = "date",
      component_pathogen_timeseries = c("alpha")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent component pathogen column
  expect_error(
    multiple(
      data = sarscov2,
      case_timeseries = "cases",
      time = "date",
      component_pathogen_timeseries = c("alpha", "nonexistent")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-numeric component pathogen column
  test_data <- sarscov2
  test_data$text_col <- as.character(test_data$alpha)
  expect_error(
    multiple(
      data = test_data,
      case_timeseries = "cases",
      time = "date",
      component_pathogen_timeseries = c("text_col")
    ),
    "numeric"
  )
})

# Tests for subtyped() function
test_that("subtyped() creates correct structure", {
  check_package_data()

  result <- subtyped(
    data = influenza,
    case_timeseries = "ili",
    time = "week",
    parent_timeseries = "inf_A",
    subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
    other_pathogen_timeseries = c("inf_B", "other")
  )

  # Test basic structure
  expect_s3_class(result, "EpiStrainDynamics.pathogen_structure")
  expect_equal(result$pathogen_structure, "subtyped")
  expect_equal(result$pathogen_names, c("inf_H3N2", "inf_H1N1", "inf_B", "other"))

  # Test data structure
  expect_equal(length(result$data$case_timeseries), nrow(influenza))
  expect_equal(nrow(result$data$component_pathogens), 3) # inf_A + inf_B + other
  expect_equal(ncol(result$data$component_pathogens), nrow(influenza))
  expect_equal(nrow(result$data$subtyped), 2) # H3N2 + H1N1
  expect_equal(ncol(result$data$subtyped), nrow(influenza))
})

test_that("subtyped() correctly orders and transposes pathogen data", {
  check_package_data()

  result <- subtyped(
    data = influenza,
    case_timeseries = "ili",
    time = "week",
    parent_timeseries = "inf_A",
    subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
    other_pathogen_timeseries = c("inf_B", "other")
  )

  # Test pathogen name ordering
  expect_equal(result$pathogen_names, c("inf_H3N2", "inf_H1N1", "inf_B", "other"))

  # Test matrix transposition for component_pathogens
  expect_equal(result$data$component_pathogens[1, ], influenza$inf_A)
  expect_equal(result$data$component_pathogens[2, ], influenza$inf_B)

  # Test matrix transposition for subtyped
  expect_equal(result$data$subtyped[1, ], influenza$inf_H3N2)
  expect_equal(result$data$subtyped[2, ], influenza$inf_H1N1)
})

test_that("subtyped() validates column existence and types", {
  check_package_data()

  # Test non-existent time column
  expect_error(
    subtyped(
      data = influenza,
      case_timeseries = "ili",
      time = "nonexistent",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent case column
  expect_error(
    subtyped(
      data = influenza,
      case_timeseries = "nonexistent",
      time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent parent column
  expect_error(
    subtyped(
      data = influenza,
      case_timeseries = "ili",
      time = "week",
      parent_timeseries = "nonexistent",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent subtyped column
  expect_error(
    subtyped(
      data = influenza,
      case_timeseries = "ili",
      time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("nonexistent"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-existent other_pathogen column
  expect_error(
    subtyped(
      data = influenza,
      case_timeseries = "ili",
      time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("nonexistent")
    ),
    "Column .nonexistent. not found in data"
  )

  # Test non-numeric columns
  test_data <- influenza
  test_data$text_col <- as.character(test_data$inf_A)

  expect_error(
    subtyped(
      data = test_data,
      case_timeseries = "text_col",
      time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "numeric"
  )

  expect_error(
    subtyped(
      data = test_data,
      case_timeseries = "ili",
      time = "week",
      parent_timeseries = "text_col",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "numeric"
  )
})

# Tests for missing data handling
test_that("all pathogen structure functions handle missing data", {
  check_package_data()

  # Create data with missing values
  test_sarscov2 <- sarscov2
  test_sarscov2$cases[5] <- NA

  test_influenza <- influenza
  test_influenza$ili[10] <- NA

  # Test single() with missing data
  expect_error(
    single(data = test_sarscov2, case_timeseries = "cases", time = "date"),
    "missing|NA"
  )

  # Test multiple() with missing data in case_timeseries
  expect_error(
    multiple(
      data = test_sarscov2,
      case_timeseries = "cases",
      time = "date",
      component_pathogen_timeseries = c("alpha")
    ),
    "missing|NA"
  )

  # Test multiple() with missing data in component pathogen
  test_sarscov2_2 <- sarscov2
  test_sarscov2_2$alpha[5] <- NA
  expect_error(
    multiple(
      data = test_sarscov2_2,
      case_timeseries = "cases",
      time = "date",
      component_pathogen_timeseries = c("alpha")
    ),
    "missing|NA"
  )

  # Test subtyped() with missing data
  expect_error(
    subtyped(
      data = test_influenza,
      case_timeseries = "ili",
      time = "week",
      parent_timeseries = "inf_A",
      subtyped_timeseries = c("inf_H3N2"),
      other_pathogen_timeseries = c("inf_B")
    ),
    "missing|NA"
  )
})

# Tests for data structure integrity
test_that("all pathogen structure functions return required components", {
  check_package_data()

  # Test single() components
  single_result <- single(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date"
  )

  expect_true(all(c(
    "pathogen_structure", "pathogen_names", "validated_tsbl",
    "data"
  ) %in% names(single_result)))
  expect_true(all(c("case_timeseries") %in% names(single_result$data)))
  expect_s3_class(single_result, "EpiStrainDynamics.pathogen_structure")

  # Test multiple() components
  multiple_result <- multiple(
    data = sarscov2,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("alpha", "delta")
  )

  expect_true(all(c(
    "pathogen_structure", "pathogen_names", "validated_tsbl",
    "data"
  ) %in% names(multiple_result)))
  expect_true(all(c("case_timeseries", "component_pathogens") %in%
    names(multiple_result$data)))
  expect_s3_class(multiple_result, "EpiStrainDynamics.pathogen_structure")

  # Test subtyped() components
  subtyped_result <- subtyped(
    data = influenza,
    case_timeseries = "ili",
    time = "week",
    parent_timeseries = "inf_A",
    subtyped_timeseries = c("inf_H3N2", "inf_H1N1"),
    other_pathogen_timeseries = c("inf_B")
  )

  expect_true(all(c(
    "pathogen_structure", "pathogen_names", "validated_tsbl",
    "data"
  ) %in% names(subtyped_result)))
  expect_true(all(c(
    "case_timeseries", "component_pathogens",
    "subtyped"
  ) %in% names(subtyped_result$data)))
  expect_s3_class(subtyped_result, "EpiStrainDynamics.pathogen_structure")
})

# Tests for edge cases
test_that("pathogen structure functions handle edge cases correctly", {
  check_package_data()

  # Test with minimal valid data (small dataset)
  small_data <- data.frame(
    date = as.Date("2020-01-01") + 0:4,
    cases = c(10, 15, 20, 25, 30),
    variant_a = c(5, 8, 10, 12, 15),
    variant_b = c(5, 7, 10, 13, 15)
  )

  result_single <- single(
    data = small_data,
    case_timeseries = "cases",
    time = "date"
  )
  expect_equal(length(result_single$data$case_timeseries), 5)

  result_multiple <- multiple(
    data = small_data,
    case_timeseries = "cases",
    time = "date",
    component_pathogen_timeseries = c("variant_a", "variant_b")
  )
  expect_equal(ncol(result_multiple$data$component_pathogens), 5)
  expect_equal(nrow(result_multiple$data$component_pathogens), 2)
})
