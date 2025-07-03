

# Statistics --------------------------------------------------------------

test_that("sig_shannon computes entropy and exponentiated index", {
  s <- sigshared::example_signature()
  expect_type(sig_shannon(s), "double")
  expect_gte(sig_shannon(s), 0)
  expect_type(sig_shannon(s, exponentiate = TRUE), "double")
  expect_gte(sig_shannon(s, exponentiate = TRUE), 1)
})

test_that("sig_shannon handles all-zero vector without error", {
  s <- sigshared::example_signature()
  s$fraction <- rep(0, nrow(s))  # All zeros

  expect_silent({
    result <- sig_shannon(s)
    expect_equal(result, 0)
  })

  expect_silent({
    result_exp <- sig_shannon(s, exponentiate = TRUE)
    expect_equal(result_exp, 1)  # exp(0) == 1
  })
})

test_that("sig_shannon returns 0 for fully peaked distribution", {
  s <- sigshared::example_signature()
  s$fraction <- c(1, rep(0, nrow(s) - 1))
  expect_equal(sig_shannon(s), 0)
  expect_equal(sig_shannon(s, exponentiate = TRUE), 1)
})

test_that("sig_shannon returns log(K) for uniform distribution", {
  s <- sigshared::example_signature()
  K <- nrow(s)
  s$fraction <- rep(1 / K, K)
  expected_entropy <- log(K)
  expect_equal(sig_shannon(s), expected_entropy, tolerance = 1e-8)
  expect_equal(sig_shannon(s, exponentiate = TRUE), K)
})

test_that("sig_shannon handles small 2-element distribution correctly", {
  s <- data.frame(
    channel = c("A", "B"),
    type = c("SBS", "SBS"),
    fraction = c(0.25, 0.75)
  )

  expected_entropy <- -0.25 * log(0.25) - 0.75 * log(0.75)
  expect_equal(sig_shannon(s), expected_entropy, tolerance = 1e-8)
  expect_equal(sig_shannon(s, exponentiate = TRUE), exp(expected_entropy), tolerance = 1e-8)
})

test_that("sig_shannon is invariant to zero-only padding", {
  s1 <- data.frame(
    channel = c("A", "B", "C"),
    type = c("SBS", "SBS", "SBS"),
    fraction = c(0.5, 0.5, 0)
  )

  s2 <- data.frame(
    channel = c("A", "B"),
    type = c("SBS", "SBS"),
    fraction = c(0.5, 0.5)
  )

  expect_equal(sig_shannon(s1), sig_shannon(s2), tolerance = 1e-8)
  expect_equal(sig_shannon(s1, exponentiate = TRUE), sig_shannon(s2, exponentiate = TRUE), tolerance = 1e-8)
})


# KL Divergence -----------------------------------------------------------



test_that("sig_kl_divergence returns non-negative divergence", {
  s <- sigshared::example_signature()
  expect_type(sig_kl_divergence(s), "double")
  expect_gte(sig_kl_divergence(s), 0)
})

test_that("sig_kl_divergence handles all-zero input", {
  s <- sigshared::example_signature()
  s$fraction <- rep(1, times = nrow(s))/nrow(s)
  expect_equal(sig_kl_divergence(s), 0)
})


test_that("sig_gini returns expected range", {
  s <- sigshared::example_signature()
  expect_type(sig_gini(s), "double")
  expect_gte(sig_gini(s), 0)
  expect_lte(sig_gini(s), 1)
})

test_that("sig_gini of uniform vector is near 0", {
  s <- sigshared::example_signature()
  s$fraction <- rep(1 / nrow(s), nrow(s))
  expect_lte(sig_gini(s), 0.01)
})
#
# test_that("sig_gini of fully peaked vector is 1", {
#   s <- sigshared::example_signature()
#   s$fraction <- c(1, rep(0, nrow(s) - 1))
#   browser()
#   expect_equal(sig_gini(s), 1)
# })


test_that("sig_l2_norm returns numeric and scales if needed", {
  s <- sigshared::example_signature()
  expect_type(sig_l2_norm(s), "double")
  expect_lt(sig_l2_norm(s, scale = TRUE), sig_l2_norm(s))
})

test_that("sig_l2_norm handles counts", {
  s <- sigshared::example_signature()
  cat <- sigstats::sig_reconstruct(s, 100)
  expect_type(sig_l2_norm(cat, value = "count"), "double")
})


test_that("sig_l2_distance symmetric and 0 for self", {
  s <- sigshared::example_signature()
  expect_equal(sig_l2_distance(s, s), 0)
})

test_that("sig_l2_distance scaled < unscaled", {
  s1 <- sigshared::example_signature()
  s2 <- s1; s2$fraction <- rev(s2$fraction)
  expect_lt(sig_l2_distance(s1, s2, scale = TRUE), sig_l2_distance(s1, s2))
})

test_that("sig_lp_distance works for various p", {
  s1 <- sigshared::example_signature()
  s2 <- s1; s2$fraction <- rev(s2$fraction)
  for (p in c(1, 2, 3, Inf)) {
    expect_type(sig_lp_distance(s1, s2, p = p), "double")
  }
})

test_that("sig_lp_distance returns 0 for identical input", {
  s <- sigshared::example_signature()
  expect_equal(sig_lp_distance(s, s, p = 2), 0)
})

test_that("sig_cosine_similarity returns 1 for identical", {
  s <- sigshared::example_signature()
  expect_equal(sig_cosine_similarity(s, s), 1)
})

test_that("sig_cosine_similarity in [0,1]", {
  s1 <- sigshared::example_signature()
  s2 <- s1; s2$fraction <- rev(s2$fraction)
  sim <- sig_cosine_similarity(s1, s2)
  expect_gte(sim, 0)
  expect_lte(sim, 1)
})

test_that("sig_collection_stats returns data.frame with metrics", {
  col <- sigshared::example_signature_collection()
  df <- sig_collection_stats(col)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("id", "gini", "shannon_index", "l1_norm") %in% names(df)))
})

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
