test_that("sig_combine works", {
  signatures = sigshared::example_valid_signature_collection()
  expect_error(sig_combine(signatures, model = c('sig1' = 0.1, 'sig2' = 0.3)), NA)
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

  # Expect cosine similarity of 1
  expect_equal(sig_cosine_similarity(sig1, sig2), expected = 1)
  expect_equal(sig_cosine_similarity(sig1, sig3), expected = 1)


  #Expect non 1
  expect_equal(sig_cosine_similarity(sig1, sig4), expected = 0.97463185)
  expect_equal(sig_cosine_similarity(sig1, sig5), expected = 0.53625854)

  # Fails if signatures have different channels/types
  expect_error(sig_cosine_similarity(sig1, sig6), regexp = "different")
  expect_error(sig_cosine_similarity(sig1, sig7), regexp = "different")
  expect_error(sig_cosine_similarity(sig1, sig8), regexp = "different")
})

