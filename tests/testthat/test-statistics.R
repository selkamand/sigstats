

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



# Gini --------------------------------------------------------------------


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

# L2 norm --------------------------------------------------------------------

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

# Lp distance --------------------------------------------------------------------
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

# Cosine distance --------------------------------------------------------------------
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


# Collection Stats --------------------------------------------------------
test_that("sig_collection_stats returns data.frame with metrics", {
  col <- sigshared::example_signature_collection()
  df <- sig_collection_stats(col)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("id", "gini", "shannon_index", "l1_norm") %in% names(df)))
})

test_that("sig_collection_stats l0_norm matches expected non-zero count", {
  sigs <- sigshared::example_signature_collection()
  stats <- sig_collection_stats(sigs)

  for (id in stats$id) {
    sig <- sigs[[id]]
    expected_l0 <- sum(sig$fraction != 0)
    expect_equal(stats[stats$id == id, "l0_norm"], expected_l0)
  }
})


# Compute Pairwise Metrics ------------------------------------------------
test_that("sig_collection_pairwise_stats returns correct format and values", {
  col <- sigshared::example_signature_collection()

  # Cosine similarity - data.frame output
  df <- sig_collection_pairwise_stats(col, metric = "cosine_similarity", format = "data.frame")
  expect_s3_class(df, "data.frame")
  expect_named(df, c("S1", "S2", "cosine_similarity"), ignore.order = TRUE)
  expect_true(all(df$cosine_similarity >= 0 & df$cosine_similarity <= 1))

  # Cosine similarity - matrix output
  m <- sig_collection_pairwise_stats(col, metric = "cosine_similarity", format = "matrix")
  expect_true(is.matrix(m))
  expect_equal(rownames(m), colnames(m))
  expect_true(all(diag(m) == 1))
  expect_true(all(m >= 0 & m <= 1, na.rm = TRUE))
  expect_equal(m, t(m))  # symmetry

  # L2 distance - data.frame
  df_l2 <- sig_collection_pairwise_stats(col, metric = "L2", format = "data.frame")
  expect_s3_class(df_l2, "data.frame")
  expect_named(df_l2, c("S1", "S2", "L2"))
  expect_true(all(df_l2$L2 >= 0))

  # L2 distance - matrix
  m_l2 <- sig_collection_pairwise_stats(col, metric = "L2", format = "matrix")
  expect_equal(rownames(m_l2), colnames(m_l2))
  expect_true(all(diag(m_l2) == 0))
  expect_equal(m_l2, t(m_l2))

  # L1 distance
  df_l1 <- sig_collection_pairwise_stats(col, metric = "L1", format = "data.frame")
  expect_named(df_l1, c("S1", "S2", "L1"))
  expect_true(all(df_l1$L1 >= 0))

  # Error for unsupported metric
  expect_error(sig_collection_pairwise_stats(col, metric = "unsupported"), "must be one of")

  # Error for bad input type
  expect_error(sig_collection_pairwise_stats(list(x = 1:10)), "not a valid signature collection")

  # Matrix input also works
  mx <- sigshared::sig_collection_reformat_list_to_matrix(col, values = "fraction")
  df_mx <- sig_collection_pairwise_stats(mx, metric = "cosine_similarity", format = "data.frame")
  expect_s3_class(df_mx, "data.frame")
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


# Model Correctness statistics --------------------------------------------

test_that("sig_model_correctness returns all expected metrics and correct types", {
  observed <- c(S1 = 0.5, S2 = 0.3, S3 = 0.2)
  truth <- c(S1 = 0.6, S2 = 0.4, S3 = 0.0)
  metrics <- sig_model_correctness(observed, truth)

  # Expect a named list with all metrics
  expect_type(metrics, "list")
  expect_true(all(c(
    "fitting_error", "RMSE", "n_false_positives", "n_true_positives",
    "n_false_negatives", "n_true_negatives", "total_false_positive_contributions",
    "precision", "recall", "specificity", "mathews_correlation_coeff",
    "f1", "balanced_accuracy"
  ) %in% names(metrics)))

  # All numeric except possibly MCC (can be NaN/NA)
  for (nm in setdiff(names(metrics), character(0))) {
    expect_true(is.numeric(metrics[[nm]]) || is.na(metrics[[nm]]))
  }
})

test_that("sig_model_correctness matches expected values for simple case", {
  observed <- c(A = 0.5, B = 0.5, C = 0)
  truth <- c(A = 1, B = 0, C = 0)
  m <- sig_model_correctness(observed, truth)
  expect_equal(m$fitting_error, 0.5)
  expect_equal(m$n_false_positives, 1)
  expect_equal(m$n_false_negatives, 0)
  expect_equal(m$n_true_positives, 1)
  expect_equal(m$n_true_negatives, 1)
  expect_equal(m$total_false_positive_contributions, 0.5)
  expect_equal(m$precision, 0.5)
  expect_equal(m$recall, 1)
  expect_equal(m$specificity, 1 / 2)
  expect_equal(m$mathews_correlation_coeff, 0.5)
})

test_that("sig_model_correctness handles perfect fit", {
  observed <- c(S1 = 0.6, S2 = 0.4, S3 = 0)
  truth <- c(S1 = 0.6, S2 = 0.4, S3 = 0)
  m <- sig_model_correctness(observed, truth)
  expect_equal(m$fitting_error, 0)
  expect_equal(m$RMSE, 0)
  expect_equal(m$n_false_positives, 0)
  expect_equal(m$n_false_negatives, 0)
  expect_equal(m$n_true_positives, 2)
  expect_equal(m$n_true_negatives, 1)
  expect_equal(m$total_false_positive_contributions, 0)
  expect_equal(m$precision, 1)
  expect_equal(m$recall, 1)
  expect_equal(m$specificity, 1)
  expect_equal(m$f1, 1)
})

test_that("sig_model_correctness handles all-zero model", {
  observed <- c(S1 = 0, S2 = 0, S3 = 0)
  truth <- c(S1 = 0.6, S2 = 0.4, S3 = 0)
  m <- sig_model_correctness(observed, truth)
  expect_equal(m$n_false_positives, 0)
  expect_equal(m$n_true_positives, 0)
  expect_equal(m$n_false_negatives, 2)
  expect_equal(m$n_true_negatives, 1)
  expect_true(is.nan(m$precision) || is.na(m$precision)) # TP=0, FP=0
  expect_equal(m$recall, 0)
  expect_equal(m$specificity, 1)
})

test_that("sig_model_correctness detects all false positives", {
  observed <- c(S1 = 0.1, S2 = 0.2, S3 = 0.7)
  truth <- c(S1 = 0, S2 = 0, S3 = 0)
  m <- sig_model_correctness(observed, truth)
  expect_equal(m$n_false_positives, 3)
  expect_equal(m$n_true_positives, 0)
  expect_equal(m$n_false_negatives, 0)
  expect_equal(m$n_true_negatives, 0)
  expect_equal(m$total_false_positive_contributions, 1)
  expect_equal(m$recall, NaN) # 0/0
  expect_equal(m$precision, 0)
  expect_equal(m$specificity, 0)
})

test_that("sig_model_correctness handles mixed ordering and missing names", {
  observed <- c(S3 = 0, S2 = 0.4, S1 = 0.6)
  truth <- c(S1 = 0.6, S2 = 0.4, S3 = 0)
  m1 <- sig_model_correctness(observed, truth)
  m2 <- sig_model_correctness(truth, observed)
  # Should be symmetric
  expect_equal(m1$fitting_error, m2$fitting_error)
  expect_equal(m1$n_true_positives, m2$n_true_positives)
})


test_that("sig_model_correctness works for unnamed vectors with validate=FALSE", {
  observed <- c(0.6, 0.4, 0)
  truth <- c(0.6, 0.4, 0)
  expect_silent(sig_model_correctness(observed, truth, validate = FALSE))
})

test_that("sig_model_correctness can handle one-hot signature (perfect/peaked)", {
  observed <- c(S1 = 1, S2 = 0, S3 = 0)
  truth <- c(S1 = 1, S2 = 0, S3 = 0)
  m <- sig_model_correctness(observed, truth)
  expect_equal(m$fitting_error, 0)
  expect_equal(m$n_true_positives, 1)
  expect_equal(m$n_false_positives, 0)
  expect_equal(m$n_false_negatives, 0)
  expect_equal(m$n_true_negatives, 2)
  expect_equal(m$precision, 1)
  expect_equal(m$recall, 1)
})

test_that("sig_model_correctness returns NA/NaN for undefined metrics", {
  # All predictions and truth are zero
  observed <- c(S1 = 0, S2 = 0)
  truth <- c(S1 = 0, S2 = 0)
  m <- sig_model_correctness(observed, truth)
  expect_true(is.nan(m$precision) || is.na(m$precision))
  expect_true(is.nan(m$recall) || is.na(m$recall))
  expect_true(is.nan(m$f1) || is.na(m$f1))
  expect_equal(m$specificity, 1)
})

test_that("sig_model_correctness is robust to floating point near-zero", {
  observed <- c(S1 = 1e-16, S2 = 1 - 1e-16)
  truth <- c(S1 = 0, S2 = 1)
  m <- sig_model_correctness(observed, truth)
  # Should treat S1 as present (since >0)
  expect_equal(m$n_false_positives, 1)
  expect_equal(m$n_true_positives, 1)
})

test_that("sig_model_correctness handles high-dimensional input", {
  n <- 100
  observed <- rep(0, n); names(observed) <- paste0("S", seq_len(n))
  truth <- observed; truth[1:5] <- 1 / 5
  m <- sig_model_correctness(observed, truth)
  expect_equal(m$n_false_negatives, 5)
  expect_equal(m$n_true_negatives, n - 5)
})

test_that("sig_model_correctness matches example in documentation", {
  observed <- c(SBS1 = 0.7, SBS5 = 0.3, SBS18 = 0)
  truth <- c(SBS1 = 0.6, SBS5 = 0.4, SBS18 = 0)
  m <- sig_model_correctness(observed, truth)
  expect_type(m, "list")
  expect_equal(m$n_true_positives, 2)
  expect_equal(m$n_true_negatives, 1)
  expect_equal(m$n_false_positives, 0)
  expect_equal(m$n_false_negatives, 0)
})


