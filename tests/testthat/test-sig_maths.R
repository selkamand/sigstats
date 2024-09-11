test_that("sig_combine works", {
  signatures <-  sigshared::example_signature_collection()
  model <- c('sig1' = 0.1, 'sig2' = 0.3)
  expect_error(sig_combine(signatures, model = model, format="combined"), NA)


  # Check format = signature outputs a valid signature object
  expect_error(
    sigshared::assert_signature(sig_combine(signatures, model = model, format = "signature"), must_sum_to_one = FALSE),
    NA
  )

  # Check Reconstruction Works with signatures where fraction sums to < 1
  combined_signature <- sig_combine(signatures, model = model, format = "signature")
  expect_error(sig_reconstruct(combined_signature, n = 100), regexp = NA)

  # Check combination works when model is 0-length numeric vector
  expect_error(sig_combine(signatures, model = numeric(0), format = "signature"), NA)
  empty_model_result <- sig_combine(signatures, model = numeric(0), format = "signature")
  expect_equal(nrow(empty_model_result), 3)
  expect_equal(empty_model_result[["signature"]], signatures[[1]][["signature"]])
  expect_equal(empty_model_result[["type"]], signatures[[1]][["type"]])
  expect_equal(empty_model_result[["fraction"]], c(0, 0, 0))

  # Fails if model is 0-length numeric vector and format=combined
  expect_error(sig_combine(signatures, model = numeric(0), format = "combined"), "sensible")

  # Check combination works when model is NULL
  expect_error(sig_combine(signatures, model = NULL, format = "signature"), NA)
  empty_model_result <- sig_combine(signatures, model = NULL, format = "signature")
  expect_equal(nrow(empty_model_result), 3)
  expect_equal(empty_model_result[["signature"]], signatures[[1]][["signature"]])
  expect_equal(empty_model_result[["type"]], signatures[[1]][["type"]])
  expect_equal(empty_model_result[["fraction"]], c(0, 0, 0))

  # Fails if model is 0-length numeric vector and format=combined
  expect_error(sig_combine(signatures, model = NULL, format = "combined"), "sensible")

  # Has some tolerance in case signature models add up to 1.000000000000000444089
  expect_error(sig_combine(signatures, model = c(sig1 = 1.000000000000000444089), format = "signature"), NA)


})

test_that("sig_cosine_similarity works", {

  # Function that just forces fraction to sum to 1
  normalize_fractions <- function(sig) {
    sig$fraction <- sig$fraction / sum(sig$fraction)
    return(sig)
  }

  sig1 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>T'),
    channel =  c('C>A', 'C>G', 'C>T'),
    fraction = c(1, 2, 3)
  ))

  sig2 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>T'),
    channel =  c('C>A', 'C>G', 'C>T'),
    fraction = c(1, 2, 3)
  ))

  sig3 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>T'),
    channel =  c('C>A', 'C>G', 'C>T'),
    fraction = c(100, 200, 300)
  ))

  sig4 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>T'),
    channel =  c('C>A', 'C>G', 'C>T'),
    fraction = c(4, 5, 6)
  ))

  sig5 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>T'),
    channel =  c('C>A', 'C>G', 'C>T'),
    fraction = c(1, 2000, 4)
  ))

  sig6 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>ADWLKAL'),
    channel =  c('C>A', 'C>G', 'C>T'),
    fraction = c(1, 2, 3)
  ))

  sig7 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G', 'C>T'),
    channel =  c('C>A', 'C>G', 'C>P'),
    fraction = c(1, 2, 3)
  ))

  sig8 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>G'),
    channel =  c('C>A', 'C>G'),
    fraction = c(1, 2)
  ))

  sig9 <- normalize_fractions(data.frame(
    type = c('C>G', 'C>A', 'C>T'),
    channel =  c('C>G', 'C>A', 'C>T'),
    fraction = c(2, 1, 3)
  ))

  sig10 <- normalize_fractions(data.frame(
    type = c('C>A', 'C>T', 'C>G'),
    channel =  c('C>A', 'C>T', 'C>G'),
    fraction = c(1, 3, 2)
  ))

  # Expect cosine similarity of 1
  expect_equal(sig_cosine_similarity(sig1, sig2), expected = 1)
  expect_equal(sig_cosine_similarity(sig1, sig3), expected = 1)


  #Expect non 1
  expect_equal(sig_cosine_similarity(sig1, sig4), expected = 0.97463185)
  expect_equal(sig_cosine_similarity(sig1, sig5), expected = 0.53625854)

  #Expect 0 when empty signatures are used
  expect_equal(sig_cosine_similarity(sigshared::example_signature_empty(), sigshared::example_signature_empty()), expected = 0)
  expect_equal(sig_cosine_similarity(sigshared::example_signature(), sigshared::example_signature_empty()), expected = 0)

  # Works when comparing sigs to catalogues
  expect_equal(sig_cosine_similarity(
    sigshared::example_signature(),
    sigshared::example_catalogue()),
    expected = 0.84672372
  )

  # Fails if signatures have different channels/types
  expect_error(sig_cosine_similarity(sig1, sig6), regexp = "different")
  expect_error(sig_cosine_similarity(sig1, sig7), regexp = "different")
  expect_error(sig_cosine_similarity(sig1, sig8), regexp = "different")

  # Works even if signatures have different sorting of channels
  expect_no_error(sig_cosine_similarity(sig1, sig9))
  expect_equal(sig_cosine_similarity(sig1, sig9), expected = 1)
  expect_equal(sig_cosine_similarity(sig9, sig10), expected = 1)
})

test_that("sig_reconstruct works", {
  signatures <-  sigshared::example_signature_collection()
  model <- c('sig1' = 0.1, 'sig2' = 0.3)
  model_sig <- sig_combine(signatures, model = model, format = "signature")
  reconstructed_catalogue <- sig_reconstruct(model_sig, n = 100)

  # Check format = signature outputs a valid signature object
  expect_error(
    sigshared::assert_catalogue(reconstructed_catalogue, must_sum_to_one = TRUE),
    NA
  )

  # Test it works with empty signatures (fractions are all 0
  sig_empty <- sigshared::example_signature()
  sig_empty[["fraction"]] <- 0

  expect_error(sig_reconstruct(signature = sig_empty, n = 100), NA)
  sig_empty_cat = sig_reconstruct(signature = sig_empty, n = 100)
  expect_equal(sig_empty_cat[["fraction"]], c(0, 0, 0))


})
