`%||%` <- function(lhs, rhs){
 if(is.null(lhs))
   return(rhs)
  else
    return(lhs)
}



# Statistics --------------------------------------------------------------
#' Calculate the Shannon Diversity Index for a Signature
#'
#' Computes the Shannon diversity index for a `sigverse` signature. This metric
#' quantifies the entropy or uncertainty associated with the distribution of mutation
#' contexts in a signature, based on the relative `fraction` of mutations in each context.
#'
#' By default, the function returns the Shannon index as an entropy value. If `exponentiate = TRUE`,
#' the function returns the **exponentiated Shannon index**, also known as the **effective number
#' of contexts** (or Hill number of order 1). This makes interpretation more intuitive:
#'
#' - A signature concentrated entirely in a single context has an exponentiated index of 1.
#' - A perfectly uniform signature (equal weight across all 96 SBS contexts) has an exponentiated index of 96.
#'
#' In biological terms, the exponentiated Shannon index answers:
#' _"How many equally frequent mutation contexts would give this level of diversity?"_
#'
#' @param signature A `sigverse` signature data.frame. See [sigshared::example_signature()].
#' @param exponentiate Logical. If `TRUE`, returns the exponentiated Shannon index (effective number of contexts).
#'
#' @return A numeric value: either the Shannon index (entropy) or the exponentiated index (effective diversity).
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#' sbs3 <- signatures[["SBS3"]]
#' sbs48 <- signatures[["SBS48"]]
#'
#' # Shannon entropy
#' sig_shannon(sbs3)
#'
#' # Exponentiated Shannon index (effective # of active contexts)
#' sig_shannon(sbs3, exponentiate = TRUE)
#'
#' # Compare with a highly focused signature
#' sig_shannon(sbs48, exponentiate = TRUE)
#'
#' @export
sig_shannon <- function(signature, exponentiate = FALSE){
  sigshared::assert_signature(signature)
  shannon_index <- compute_shannon_index(signature[["fraction"]])

  if(exponentiate){
   return(exp(shannon_index))
  }

  return(shannon_index)
}

compute_shannon_index <- function(probabilities, exponentiate = FALSE){
  shannon_index <- -sum(probabilities * log(probabilities))


  if(exponentiate){
    return(exp(shannon_index))
  }

  return(shannon_index)
}


#' Compute KL Divergence from a Uniform Signature
#'
#' Calculates the Kullback–Leibler (KL) divergence between a mutational signature
#' and a uniform distribution. KL divergence quantifies how much the observed
#' mutation context distribution (the signature) deviates from an equal-weight,
#' flat profile across all contexts.
#'
#' A value of 0 indicates a perfectly uniform signature. Higher values indicate
#' more peaked or biased signatures. KL divergence is commonly used as a measure
#' of "non-uniformity" or "distinctiveness" of mutation profiles.
#'
#' A small `pseudocount` is added to avoid taking the log of zero when any context
#' has zero weight.
#'
#' @param signature A `sigverse` signature data.frame. Must contain a `fraction` column.
#' @param base The logarithmic base to use (default is natural log, `exp(1)`). Use 2 for bits.
#' @param pseudocount A small positive number added to each term to prevent log(0).
#'
#' @return A single numeric value representing the KL divergence from uniform.
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#' sbs3 <- signatures[["SBS3"]]
#'
#' # Compute KL divergence (how far is SBS3 from flat?)
#' sig_kl_divergence(sbs3)
#'
#' # Use base 2 (bits)
#' sig_kl_divergence(sbs3, base = 2)
#'
#' @seealso [sig_shannon()], [sig_gini()]
#' @export
sig_kl_divergence <- function(signature, base = exp(1), pseudocount = 1e-12){

  # Assertions
  sigshared::assert_signature(signature)

  observed <- signature[["fraction"]]
  expected <- rep(1,times = length(observed))/length(observed)

  # Computation
  compute_kl_divergence(
    p = observed,
    q = expected,
    pseudocount = pseudocount,
    base = base
  )
}

compute_kl_divergence <- function(p, q, pseudocount = 1e-12, base = exp(1)){
  if(all(p == 0)) return(0)
  sum(p*log(p/(q+pseudocount), base = base), na.rm = TRUE)
}


#' Compute the Gini Coefficient of a Signature or Catalogue
#'
#' Calculates the **Gini coefficient**, a measure of inequality or concentration,
#' for a `sigverse` signature or catalogue. It ranges from:
#'
#' - **0**: perfectly uniform distribution (e.g. all 96 contexts equal)
#' - **1**: total concentration in a single context
#'
#' The Gini coefficient complements entropy-based measures like Shannon index by capturing
#' the *inequality* of the distribution, rather than uncertainty or diversity.
#'
#' This function uses the **biased version** of the Gini coefficient (dividing by _n_), which is:
#' - The **standard definition** used in descriptive analysis
#' - More appropriate when the signature is a full population profile rather than a random sample
#' - Consistent with usage in mutation signature literature
#'
#' @param signature A `sigverse` signature or catalogue data.frame.
#'
#' @return A numeric value between 0 and 1 representing the Gini coefficient.
#'
#' @examples
#' library(sigstash)
#' sigs <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' sig_gini(sigs[["SBS1"]])  # moderately peaked
#' sig_gini(sigs[["SBS48"]]) # highly peaked, close to 1
#'
#'
#' @seealso [sig_shannon()], [sig_kl_divergence()], [sig_l2_norm()]
#' @export
sig_gini <- function(signature){
  sigshared::assert_signature(signature)

  vec <- signature[["fraction"]]
  compute_gini(vec)
}

compute_gini <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (all(x == 0)) return(0)  # define Gini = 0 when all entries are 0
  G <- sum((2 * seq_len(n) - n - 1) * x)
  G / (n * sum(x))
}

#' Compute the L2 Norm of a Signature or Catalogue
#'
#' Calculates the **L2 norm** (Euclidean norm) of a `sigverse` signature or catalogue.
#' This provides a quantitative measure of how concentrated or dispersed the distribution is.
#'
#' For a vector \( x \), the L2 norm is defined as:
#' \deqn{\|x\|_2 = \sqrt{\sum_i x_i^2}}
#'
#' Interpretation:
#' - A signature with a **uniform distribution** has a **lower** L2 norm.
#' - A **peaked signature** (one dominant context) has a **higher** L2 norm.
#'
#' This is an alternative to entropy-based metrics (like Shannon index), and useful
#' for quickly identifying signatures with strong focal points.
#'
#' @param signature A `sigverse` signature or catalogue data.frame.
#' @param value Character string, either `"fraction"` or `"count"`, indicating which column to compute the norm on.
#' @param scale Logical. If `TRUE`, divides the norm by the number of elements to enable easier comparison across different signature sizes.
#'
#' @return A single numeric value representing the L2 norm.
#'
#' @examples
#' library(sigstash)
#'
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Compute L2 norm on fractional signature
#' sig_l2_norm(signatures[["SBS1"]])
#'
#' # Compare with a flatter signature
#' sig_l2_norm(signatures[["SBS3"]])
#'
#' # Compute on raw counts (requires a catalogue)
#' cat1 <- sig_reconstruct(signatures[["SBS3"]], n = 100)
#' sig_l2_norm(cat1, value = "count")
#'
#' @seealso [sig_shannon()], [sig_kl_divergence()], [sig_gini()]
#' @export
sig_l2_norm <- function(signature, value = c("fraction", "count"), scale = FALSE){

  # Assertions
  value <- rlang::arg_match(value)
  requires_catalogue <- value == "count"
  sigshared::assert_signature(signature, allow_catalogues = requires_catalogue)

  # Computation
  vec <- signature[[value]]
  norm <- compute_l2_norm(vec)

  # Scale by number of channels
  if(scale) {
    norm <- norm/length(vec)
  }

  return(norm)
}


#' Compute the L2 Distance Between Two Signatures or Catalogues
#'
#' Calculates the **L2 distance** (Euclidean distance) between two `sigverse` signatures or catalogues.
#' This metric quantifies how different the distributions are in terms of their numeric values.
#'
#' For vectors \( x \) and \( y \), the L2 distance is:
#' \deqn{\|x - y\|_2 = \sqrt{\sum_i (x_i - y_i)^2}}
#'
#' A smaller value indicates more similar signatures, while larger values indicate greater dissimilarity.
#' This is often used as a simple and fast alternative to cosine similarity or KL divergence.
#'
#' @param signature1,signature2 Two `sigverse` signatures or catalogues.
#'   See [sigshared::example_signature()] or [sigshared::example_catalogue()].
#'   Must contain matching `channel` values in identical order.
#' @param value A character string: `"fraction"` for normalised signatures or `"count"` for raw catalogues.
#' @param scale Logical. If `TRUE`, divides the L2 distance by the number of elements.
#'   This can help normalise distance values across signatures with different dimensions.
#'
#' @return A single numeric value representing the L2 distance.
#'
#' @examples
#' library(sigstash)
#'
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#' s1 <- signatures[["SBS1"]]
#' s2 <- signatures[["SBS5"]]
#'
#' # Compute distance between two fractional signatures
#' sig_l2_distance(s1, s2)
#'
#' # Compare catalogue-level distance (on raw counts)
#' cat1 <- sig_reconstruct(s1, n = 100)
#' cat2 <- sig_reconstruct(s2, n = 100)
#' sig_l2_distance(cat1, cat2, value = "count")
#'
#' @export
sig_l2_distance <- function(signature1, signature2, value = c("fraction", "count"), scale = FALSE){

  # Assertions
  value <- rlang::arg_match(value)
  requires_catalogue <- value == "count"
  sigshared::assert_signature(signature1, allow_catalogues = requires_catalogue)
  sigshared::assert_signature(signature2, allow_catalogues = requires_catalogue)
  assertions::assert_identical(signature1[["channel"]], signature2[["channel"]])

  # Computation
  vec <- signature1[[value]] - signature2[[value]]

  norm <- compute_l2_norm(vec)

  # Scale by number of channels
  if(scale) {
    norm <- norm/length(vec)
  }

  return(norm)
}

compute_l2_norm <- function(vec){
  compute_norm(vec, 2)
}

compute_norm <- function(vec, p){
  (sum(vec^p))^(1/p)
}

#' Calculate Cosine Matrix Between Two Signatures
#'
#' Computes cosine similarity between each pair of signatures in a sigverse signature collection
#'
#' @param signature1,signature2 sigverse signature data.frames. See [sigshared::example_signature()].
#' @param assume_sensible_input A logical flag indicating whether to skip validation checks for the input signatures.
#' Enabling this option can significantly speed up the cosine similarity calculation by assuming that both inputs are valid
#' signature objects and that their channels are already sorted. This option should only be used when performance is critical and these assumptions can be verified upstream.
#' @return a number between 0 and 1 representing cosine similarity
#'
#' @export
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Compute cosine similarity between two signatures
#' sig_cosine_similarity(signatures[["SBS1"]], signatures[["SBS2"]])
#'
sig_cosine_similarity <- function(signature1,signature2, assume_sensible_input = FALSE){

  if(!assume_sensible_input){
    sigshared::assert_signature(signature1, must_sum_to_one = FALSE)
    sigshared::assert_signature(signature2, must_sum_to_one = FALSE)

    # Ensure signatures are sorted the same way 'type channel' match 1:1 (including order)
    sig1_type_channel_id = paste0('type:', signature1[['type']],', ', 'channel:',signature1[['channel']])
    sig2_type_channel_id = paste0('type:', signature2[['type']],', ', 'channel:',signature2[['channel']])
    identical = identical(sig1_type_channel_id,sig2_type_channel_id)

    if(!identical){
      # Types + Channels combinations are either not identical, or sorted differently
      unique_to_sig1 <- setdiff(sig1_type_channel_id, sig2_type_channel_id)
      unique_to_sig2 <- setdiff(sig2_type_channel_id, sig1_type_channel_id)

      # Check if they're just sorted differently
      set_equivalent <- setequal(sig1_type_channel_id, sig2_type_channel_id) & length(sig1_type_channel_id) == length(sig2_type_channel_id)

      # If signatures have different channels/types throw an error
      assertions::assert(set_equivalent, msg = c(
        "Cannot calculate cosine similarity between signatures/catalogues with different types or channels:",
        c(
          "*"="Unique Types/Channels in Signature1: [{unique_to_sig1}]",
          "*"="Unique Types/Channels in Signature2: [{unique_to_sig2}]"
        )

      ))

      # If signatures are just sorted differently, fix the sorting and continue
      new_order_for_sig2 <- match(sig1_type_channel_id, sig2_type_channel_id)
      signature2 <- signature2[new_order_for_sig2,]
    }

  }

  cosine <- sim_cosine(signature1[['fraction']], signature2[['fraction']])

  # If NaN (e.g.if all fractions = 0) because replace with 0
  if(is.nan(cosine)) cosine <- 0

  return(cosine)
}


sim_cosine <- function(x, y){
  as.numeric(lsa::cosine(x, y))
}


# Operations --------------------------------------------------------------


#' Add Signatures to a Combined Signature Model
#'
#' This function takes a signature collection and a model as input and adds
#' selected signatures to a combined signature model.
#' The default output is a data.frame in either sigverse signature format or the 'combined_signature_model'
#' style depending on `format`.
#'
#' @param signatures A sigverse signature collection. See ([sigshared::example_signature_collection()]) for details.
#' @param model A named numeric vector representing the contribution of each signature to the combined model.
#' The names correspond to the signatures, and the values represent their contributions.
#' The sum of the values in this vector should be less than or equal to 1.
#' @param format A character string indicating the output format.
#' If "combined", the function returns a 'combined_signature_model' data.frame where each row represents a contribution for a particular channel from a single signature (duplicate channels are not collapsed).
#' If "signature", the function returns the data in the sigverse signature format, representing a novel signature created by combining the signatures in the collection according to the ratios described by the model.
#' @param verbose enables detailed output messages (flag).
#' @return A data.frame in the 'combined_signature_model' style containing the selected signatures and their modified fractions based on the model.
#'
#' @seealso \code{\link{example_signature_collection}} \code{\link{example_model}}
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Create a model that represents a mix of SBS1 (40%) and SBS2 (60%)
#' model <- c(SBS1 = 0.4, SBS2 = 0.6)
#'
#' # Add selected signatures to the combined model
#' combined_signatures <- sig_combine(signatures, model)
#' print(combined_signatures)
#'
#'
#' @export
sig_combine <- function(signatures, model, format = c("signature", "combined"), verbose=FALSE){

  # Settings
  tolerance <- 5e-07 # Allows models to sum to anything < 1.000001

  # Null replacements
  model <- model %||% numeric(0)

  # Assertions
  sigshared::assert_signature_collection(signatures)
  sigshared::assert_model(model, signatures)
  format <- rlang::arg_match(format)

  model_signatures <- names(model)

  # Deal with the case that model vector has no length.
  # If format="signature" we can just return an empty signature (all fraction values = 0)
  # If the format="combined" (individual signature contributions preserved in dataframe) there is no
  # Sensible way to represent an model with no contributions, and so an error is thrown
  if(length(model_signatures) == 0){
    assertions::assert(
      format == "signature",
      msg = "There is no sensible way to represent an model with no signature contributions if {.arg format='{format}'}. Try setting {.arg format='signature'}"
    )
    if(verbose) { warning("model vector supplied to sig_combine is empty. Returning an empty (fraction=0) signature") }
    df_signatures = signatures[[1]]
    df_signatures <- df_signatures[c("type", "channel", "fraction")]
    df_signatures[["fraction"]] <- 0
    return(df_signatures)
  }


  # Grab relevant signatures, modify fraction to reflect the signatures total contribution, and add signature name
  ls_signatures <- lapply(model_signatures, function(signame){
    contribution = model[signame]
    orig_sig <- signatures[[signame]]
    new_sig <- orig_sig
    new_sig[['fraction_original']] <- new_sig[['fraction']]
    new_sig[['fraction']] <- orig_sig[['fraction']] * contribution
    new_sig[['signature']] <- signame
    new_sig <- new_sig[c('signature', names(new_sig)[-5])]
    return(new_sig)
  })

  # Create a data.frame
  df_signatures_combined <- do.call(rbind, ls_signatures)

  # Format the data.frame
  if(format == "signature"){
    df_signatures_collapsed <- sig_combine_collapse_to_single_signature(df_signatures_combined)
    return(df_signatures_collapsed)
  }

  # Return the data.frame
  return(df_signatures_combined)
}

#' Collapse to single signatures
#'
#' Convert the output of [sig_combine()] into a simple sigverse signature object.
#' This is useful when you want to run maths on a signature derived from a model
#' Note for visualisation of signature combination models,
#' we suggest directly using the output of [sig_combine()] in [sigvis::sig_visualise()] so that
#' when [sigvis::sig_make_interactive()] is run the original signature contributions are preserved
#'
#' @param signature_combination a dataframe produced by [sig_combine()]
#' which represents the combination of multiple signatures with known (exact) proportions - where each individual signature is kept distinct to make it easy to plot as a stackedbar
#'
#' @return a data.frame in sigverse signature format
#' @export
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Create a model that represents a mix of SBS1 (40%) and SBS2 (60%)
#' model <- c(SBS1 = 0.4, SBS2 = 0.6)
#'
#' # Add selected signatures to the combined model
#' combined_signatures <- sig_combine(signatures, model)
#'
#' # Flatten the combined_signatures dataframe that keeps separate signatures
#' signature <- sig_combine_collapse_to_single_signature(combined_signatures)
#'
sig_combine_collapse_to_single_signature <- function(signature_combination){
  stats::aggregate(data = signature_combination, fraction ~ type + channel, FUN = function(frac){ sum(frac)})
}



#' Reconstruct a catalogue from a signature
#'
#' Generates reconstructed substitution profile (catalogues) generated by proportionally combining each estimated contribution.
#'
#' @param signature a data.frame in sigverse signature format. See [sigshared::example_signature()] for details.
#' @param n total number of variants to use for reconstruction
#' @return data.frame conforming to sigverse catalogue format
#' @export
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Create a model that represents a mix of SBS1 (40%) and SBS2 (60%)
#' model <- c(SBS1 = 0.4, SBS2 = 0.6)
#'
#' # Add selected signatures to the combined model
#' combined_signatures <- sig_combine(signatures, model)
#'
#' # Flatten the combined_signatures dataframe that keeps separate signatures
#' signature <- sig_combine_collapse_to_single_signature(combined_signatures)
#'
#' # Reconstruct a perfect catalogue describing what the mutational profile of a sample
#' # with 200 mutations and the given signature model would look like
#' reconstuction <- sig_reconstruct(signature, n=200)
#'
#' # Visualise Result
#' # sig_visualise(reconstuction, class = "catalogue")
#'
sig_reconstruct <- function(signature, n){

  # assertions
  assertions::assert_number(n)
  sigshared::assert_signature(signature, must_sum_to_one = FALSE)

  # Convert signature to catalogue
  signature[['count']] <- signature[['fraction']] * n

  # Automatically fix catalogue fraction just in case signature it was generated from doesn't sum to 1
  signature[['fraction']] <- signature[['count']]/sum(signature[['count']], na.rm = TRUE)

  # Replace NaNs with 0 in case total counts are 0 (e.g. if empty signature is supplied)
  signature[['fraction']][is.na(signature[['fraction']])] <- 0

  return(signature)
}


#' Subtract Signatures/Catalogues
#'
#' Subtracts signature2 from signature1 and returns result.
#'
#' @param signature1,signature2 sigverse signature/catalogue data.frames. See [sigshared::example_signature()] or [sigshared::example_catalogue()].
#'
#' @return a data.frame representing a sigverse signature.
#' @export
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Subtract signatures
#' signatures[['SBS3']] %-% signatures[['SBS4']]
#'
#' # Identical approach using full function name
#' sig_subtract(signatures[['SBS3']], signatures[['SBS4']])
#'
sig_subtract <- function(signature1, signature2) {
  sigshared::assert_signature(signature1, must_sum_to_one = FALSE)
  sigshared::assert_signature(signature2, must_sum_to_one = FALSE)

  # Ensure matching structure and ordering
  id1 <- paste(signature1[['type']], signature1[['channel']])
  id2 <- paste(signature2[['type']], signature2[['channel']])
  assertions::assert_equal(id1, id2, msg = "To subtract signatures they must have identical type and channel orders.")

  is_catalogue <- "count" %in% colnames(signature1) && "count" %in% colnames(signature2)

  sig_result <- signature1

  if (is_catalogue) {
    sig_result[["count"]] <- signature1[["count"]] - signature2[["count"]]
    sig_result[["fraction"]] <- compute_fraction_from_count(sig_result[["count"]])
  } else {
    sig_result[["fraction"]] <- signature1[["fraction"]] - signature2[["fraction"]]
  }

  return(sig_result)
}

#' @rdname sig_subtract
#' @export
`%-%` <- sig_subtract

#' Add Two Catalogues
#'
#' Sums two sigverse-style catalogues element-wise by their `count` values.
#'
#' This operator provides a concise way to add two catalogues. For summing multiple
#' catalogues, use [sig_sum()].
#'
#' @param catalogue1,catalogue2 Two `sigverse` catalogue data.frames.
#'  See [sigshared::example_catalogue()].
#'  Must have the same `type` and `channel` rows in the same order.
#'
#' @return A single `sigverse` catalogue data.frame with summed `count`
#'   and recomputed `fraction` values. See [sigshared::example_catalogue()].
#'
#' @seealso [sig_sum()] for summing a collection of catalogues
#'
#' @export
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#' cat1 <- sig_reconstruct(signatures[['SBS3']], n = 100)
#' cat2 <- sig_reconstruct(signatures[['SBS4']], n = 100)
#'
#' # Add two catalogues
#' cat_sum <- cat1 %+% cat2
#'
#' # Add multiple catalogues via repeated %+%
#' cat3 <- sig_reconstruct(signatures[['SBS5']], n = 100)
#' total <- cat1 %+% cat2 %+% cat3
#'
#' # Addition two signatures using the full function name
#' sig_add(cat1, cat2)
#'
#' # Alternatively, use sig_sum for a list of catalogues
#' catalogues <- list(cat1 = cat1, cat2 = cat2, cat3 = cat3)
#' total2 <- sig_sum(catalogues)
sig_add <- function(catalogue1, catalogue2){
  sigshared::assert_catalogue(catalogue1, must_sum_to_one = FALSE)
  sigshared::assert_catalogue(catalogue2, must_sum_to_one = FALSE)

  # Ensure matching structure and ordering
  id1 <- paste(catalogue1[['type']], catalogue1[['channel']])
  id2 <- paste(catalogue2[['type']], catalogue2[['channel']])
  assertions::assert_equal(id1, id2, msg = "To sum signatures they must have identical type and channel orders.")

  # Catalogue
  cat_result <- catalogue1
  cat_result[["count"]] <- catalogue1[["count"]] + catalogue2[["count"]]
  cat_result[["fraction"]] <- compute_fraction_from_count(cat_result[["count"]])

  return(cat_result)
}

#' @rdname sig_add
#' @export
`%+%` <- sig_add

#' Sum a Collection Catalogues
#'
#' Sums a list of sigverse catalogues, or two individual catalogues, into a single result.
#'
#' This function is useful for aggregating catalogues across samples or replicates.
#' If you only need to add two catalogues, you may use the [`%+%`] operator instead.
#'
#' @param catalogues A named list of `sigverse` catalogue data.frames. See [sigshared::example_catalogue_collection()]
#'
#' @return A `sigverse` catalogue data.frame representing the total. See [sigshared::example_catalogue()].
#'
#' @seealso [`%+%`] for summing two catalogues
#'
#' @export
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Reconstruct catalogues for two pure samples (each with 100 mutations)
#' catalogue1 <- sig_reconstruct(signatures[['SBS3']], n = 100)
#' catalogue2 <- sig_reconstruct(signatures[['SBS4']], n = 100)
#'
#' # Sum catalogue1 and  catalogue2
#' catalogue_sum <- catalogue1 %+% catalogue2
#'
#' # Sum a collection
#' collection <- list(cat1 = catalogue1, cat2 = catalogue2)
#' collection <- sig_sum(collection)
sig_sum <- function(catalogues){
  sigshared::assert_catalogue_collection(catalogues)
  mx <- sigshared::sig_collection_reformat_list_to_matrix(catalogues, values = "count")
  sums <- rowSums(mx)
  channels <- rownames(mx)
  types <- attr(mx, "types")

  cat_result <- catalogues[[1]]
  cat_result[["count"]] <- sums
  cat_result[["fraction"]] <- compute_fraction_from_count(cat_result[["count"]])
  return(cat_result)
}

