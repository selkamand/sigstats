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
