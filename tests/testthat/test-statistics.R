
# Computing Experimental P Value ------------------------------------------

test_that("Proportion of bootstraps below threshold is calculated correctly", {
  bootstraps <- c(0.01, 0.03, 0.02, 0.05, 0.07)
  threshold <- 0.05

  expect_equal(sig_compute_experimental_p_value(bootstraps, threshold), 0.6)
})

test_that("Function returns 0 when no bootstrap values are below the threshold", {
  bootstraps <- c(0.1, 0.2, 0.3)
  threshold <- 0.05

  expect_equal(sig_compute_experimental_p_value(bootstraps, threshold), 0)
})

test_that("Function returns 1 when all bootstrap values are below the threshold", {
  bootstraps <- c(0.01, 0.02, 0.03)
  threshold <- 0.05

  expect_equal(sig_compute_experimental_p_value(bootstraps, threshold), 1)
})

test_that("Function handles empty vector correctly", {
  bootstraps <- numeric(0)
  threshold <- 0.05

  expect_true(is.nan(sig_compute_experimental_p_value(bootstraps, threshold)))
})

test_that("Function handles vector with all equal values correctly", {
  bootstraps <- rep(0.05, 5)
  threshold <- 0.05

  expect_equal(sig_compute_experimental_p_value(bootstraps, threshold), 0)
})

test_that("Function handles different threshold values correctly", {
  bootstraps <- c(0.01, 0.03, 0.02, 0.05, 0.07)

  expect_equal(sig_compute_experimental_p_value(bootstraps, 0.01), 0)
  expect_equal(sig_compute_experimental_p_value(bootstraps, 0.07), 0.8)
  expect_equal(sig_compute_experimental_p_value(bootstraps, 0.06), 0.8)
})

test_that("Function handles large input vectors correctly", {
  bootstraps <- rnorm(10000, mean = 0.05, sd = 0.01)
  threshold <- 0.05

  proportion <- sum(bootstraps < threshold) / length(bootstraps)

  expect_equal(sig_compute_experimental_p_value(bootstraps, threshold), proportion)
})



# Summarising Bootstrap Format --------------------------------------------

test_that("sig_summarise_bootstraps handles valid input", {
  # Generate example bootstraps from sigshared
  bootstraps <- sigshared::example_bootstraps()

  # Run the summarise function with a reasonable threshold
  result <- sig_summarise_bootstraps(bootstraps, threshold = 0.05)

  # Check the structure of the result
  expect_true(is.data.frame(result))
  expect_named(result, c('iqr', 'max', 'median', 'min', 'n', 'outlier_high_threshold', 'outlier_low_threshold', 'outliers', 'p_value', 'q1', 'q3', 'signatures'), ignore.order = TRUE)

  # Check if the signatures are correct
  expect_equal(sort(unique(bootstraps$signature)), sort(result$signatures))

  # Ensure p-values are numeric and within 0 and 1
  expect_true(all(result$p_value >= 0 & result$p_value <= 1))

  # Check if outliers column is correctly formatted (default list column)
  expect_true(is.list(result$outliers))
})

test_that("sig_summarise_bootstraps handles outliers as strings", {
  # Generate example bootstraps
  bootstraps <- sigshared::example_bootstraps()

  # Run the summarise function with outliers as strings
  result <- sig_summarise_bootstraps(bootstraps, threshold = 0.05, outliers_as_strings = TRUE)

  # Ensure outliers are returned as character strings, not lists
  expect_true(is.character(result$outliers))
  expect_true(all(grepl("|", result$outliers))) # Check for '|' as a separator
})

test_that("sig_summarise_bootstraps handles invalid input", {
  # Invalid bootstraps input (not a data frame)
  expect_error(sig_summarise_bootstraps(list()), "must be represented as a data.frame")

  # Missing required columns
  bad_bootstraps <- data.frame(bootstrap = c(1, 1, 2), signature = c("Signature1", "Signature2", "Signature1"))
  expect_error(sig_summarise_bootstraps(bad_bootstraps), regexp = "must contain the following", fixed=TRUE)

})

test_that("sig_summarise_bootstraps returns correct quantiles and IQR", {
  # Generate example bootstraps
  bootstraps <- data.frame(
    bootstrap = rep(c(1, 2), each = 3L),
    signature = rep(c("Signature1", "Signature2", "Signature3"), 2),
    contribution_absolute = c(300, 690, 10, 440, 500, 60),
    contribution = c(0.3, 0.69, 0.01, 0.44, 0.5, 0.06)
  )

  # Run the summarise function
  result <- sig_summarise_bootstraps(bootstraps, threshold = 0.05)

  # Check if quantiles and IQR are correct (need to mock or compute them manually)
  summary_stats <- boxplotstats::calculate_boxplot_stats_for_multiple_groups(bootstraps$contribution, ids = bootstraps$signature)
  expect_equal(result$quantiles, summary_stats$quantiles)
  expect_equal(result$iqr, summary_stats$iqr)

  # Check column names
  expect_named(result,
               c('iqr', 'max', 'median', 'min', 'n', 'outlier_high_threshold', 'outlier_low_threshold', 'outliers', 'p_value', 'q1', 'q3', 'signatures'),
               ignore.order = TRUE
             )


  expect_equal(result$p_value, c(0, 0, 0.5))
})
