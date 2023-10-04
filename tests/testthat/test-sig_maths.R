test_that("sig_combine works", {
  signatures = sigshared::example_valid_signature_collection()
  expect_error(sig_combine(signatures, model = c('sig1' = 0.1, 'sig2' = 0.3)), NA)
})
