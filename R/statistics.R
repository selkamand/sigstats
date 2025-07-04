
#  Statistics --------------------------------------------------------------

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
  sigshared::assert_signature(signature, must_sum_to_one = FALSE)
  shannon_index <- compute_shannon_index(signature[["fraction"]])

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




#' Compute the Gini Coefficient of a Signature or Catalogue
#'
#' Calculates the **Gini coefficient**, a measure of inequality or concentration,
#' for a `sigverse` signature or catalogue. It quantifies how unevenly the mutation
#' probability mass is distributed across contexts:
#'
#' - **0**: perfectly uniform distribution (e.g. all 96 contexts equal)
#' - **1**: total concentration in a single context
#'
#' The Gini coefficient complements entropy-based measures like the Shannon index by capturing
#' the *inequality* of the distribution, rather than its uncertainty or diversity.
#'
#' This function uses the **unbiased version** of the Gini coefficient, scaled by `K / (K - 1)`
#' where `K` is the number of mutation contexts. This adjustment:
#' - Ensures the Gini ranges from 0 to 1 for any number of contexts
#' - Makes the measure **comparable across signatures** with different numbers of mutation types
#' - Is appropriate for full probability distributions (as in mutation signatures)
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
#' @seealso [sig_shannon()], [sig_kl_divergence()], [sig_l2_norm()]
#' @export
sig_gini <- function(signature){
  sigshared::assert_signature(signature)

  vec <- signature[["fraction"]]
  compute_gini(vec)
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
  norm <- compute_norm(vec, p = 2)

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
  sig_lp_distance(signature1 = signature1, signature2 = signature2, p = 2, value = value, scale = scale)
}

#' Compute L<sub>p</sub> Distance Between Two Signatures
#'
#' Calculates the L<sub>p</sub> distance (also known as the Minkowski distance)
#' between two `sigverse` signatures or catalogues. This generalizes various distance
#' metrics depending on the choice of `p`.
#'
#'
#' This function is useful for flexible distance computations when comparing
#' mutational signatures or catalogues. All channels must match and be in the same order.
#'
#' By default, distances are computed using raw values. If `scale = TRUE`, the distance is
#' divided by the number of mutation contexts to allow comparisons across different signature types
#' (e.g. SBS vs DBS).
#'
#' @param signature1, signature2 Two `sigverse` signature or catalogue data.frames.
#' @param p A numeric value ≥ 0 indicating the order of the L<sub>p</sub> norm to compute.
#'   - `p = 0`: counts the number of non-zero entries (not a true norm; useful for sparsity).
#'   - `p = 1`: Manhattan distance (sum of absolute differences).
#'   - `p = 2`: Euclidean (L2) distance.
#'   - `p = ∞`: Chebyshev distance (maximum absolute difference).
#' @param value Either `"fraction"` (default) or `"count"` — determines which column is used for comparison.
#' @param scale Logical. If `TRUE`, distance is divided by the number of contexts (i.e. channels).
#'
#' @return A non-negative numeric value representing the L<sub>p</sub> distance between the two profiles.
#'
#' @examples
#' library(sigstash)
#' sigs <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' sig_lp_distance(sigs[["SBS1"]], sigs[["SBS5"]], p = 1)  # L1 (Manhattan)
#' sig_lp_distance(sigs[["SBS1"]], sigs[["SBS5"]], p = 2)  # L2 (Euclidean)
#' sig_lp_distance(sigs[["SBS1"]], sigs[["SBS5"]], p = 1, scale = TRUE)
#'
#' @seealso [sig_l2_distance()], [sig_cosine_similarity()]
#' @export
sig_lp_distance <- function(signature1, signature2, p, value = c("fraction", "count"), scale = FALSE){

  # Assertions
  value <- rlang::arg_match(value)
  requires_catalogue <- value == "count"
  sigshared::assert_signature(signature1, allow_catalogues = requires_catalogue)
  sigshared::assert_signature(signature2, allow_catalogues = requires_catalogue)
  assertions::assert_identical(signature1[["channel"]], signature2[["channel"]])

  # Computation
  distance <- compute_lp_distance(
    signature1[[value]],
    signature2[[value]],
    p = p
  )

  # Scale by number of channels
  if(scale) {
    distance <- distance/nrow(signature1)
  }

  return(distance)
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

  cosine <- compute_cosine_similarity(signature1[['fraction']], signature2[['fraction']])

  # If NaN (e.g.if all fractions = 0) because replace with 0
  if(is.nan(cosine)) cosine <- 0

  return(cosine)
}


#' Compute Summary Statistics for a Signature Collection
#'
#' Calculates a panel of descriptive statistics for each signature in a `sigverse` collection.
#' These metrics quantify different aspects of a mutational signature's shape — such as inequality,
#' diversity, concentration, and sparsity - and are useful for comparing how "flat", "focal", or
#' "distinctive" different signatures are.
#'
#' For each signature, the following metrics are reported:
#'
#' - **`gini`**: Measures inequality (0 = perfectly flat; 1 = all weight in one context).
#' - **`shannon_index`**: Entropy of the distribution (higher = more uncertain/diverse).
#' - **`shannon_index_exp`**: Effective number of active contexts (e.g., 96 = flat, 1 = peaked).
#' - **`shannon_index_exp_scaled`**: Fraction of maximum possible diversity (0–1 scale).
#' - **`kl_divergence_from_uniform`**: Divergence from a uniform (flat) distribution.
#' - **`l1_norm`**: Total absolute weight (larger = more mass in fewer contexts).
#' - **`l2_norm`**: Magnitude of the vector; emphasizes focal peaks.
#' - **`l3_norm`**: Amplifies concentration even more than L2.
#' - **`l0_norm`**: Number of non-zero contexts (also known as the L0 "norm").
#'     This is not a true mathematical norm but is commonly used as a measure of **sparsity** —
#'     how many mutation channels contribute at all. A value of 0 means the signature is completely empty;
#'     a higher value indicates more active contexts.
#' - **`*_scaled` variants**: Norms divided by number of contexts to allow cross-signature comparisons.
#' - **`max_channel_fraction`**: Highest single-context weight (equivalent to the infinity norm).
#'
#' This function is optimised for speed (is faster than computing each norm independently for each signature)
#' and returns a data.frame with one row per signature and columns for each computed metric.
#'
#' @param signatures A `sigverse` signature collection (named list of signature data.frames).
#'
#' @return A `data.frame` with one row per signature and columns for each computed metric.
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Compute statistics for all signatures
#' stats <- sig_collection_stats(signatures)
#' head(stats)
#'
#' # Examine metrics for a single signature
#' stats[stats$id == "SBS1", ]
#'
#' @seealso [sig_shannon()], [sig_gini()], [sig_kl_divergence()], [sig_l2_norm()]
#' @export
sig_collection_stats <- function(signatures){

  # Assertions
  sigshared::assert_signature_collection(signatures)

  # Convert to matrix (makes computation faster)
  mx_signatures <- sigshared::sig_collection_reformat_list_to_matrix(signatures, values = "fraction")

  # Compute uniform distribution (used later for kl divergence calculation)
  unif = rep(1, times = nrow(mx_signatures))/nrow(mx_signatures)

  # Scale norms based on number of channels if scale = TRUE
  n_channels <- nrow(mx_signatures)

  # Metrics
  ls_metrics <- apply(
    X = mx_signatures,
    MARGIN = 2,
    function(x){
      c(
        gini = compute_gini(x),
        shannon_index = compute_shannon_index(x, exponentiate = FALSE),
        shannon_index_exp = compute_shannon_index(x, exponentiate = TRUE),
        kl_divergence_from_uniform = compute_kl_divergence(p = x, q = unif),
        l3_norm = compute_norm(x, p = 3),
        l2_norm = compute_norm(x, p = 2), # Euclidean Distance
        l1_norm = compute_norm(x, p = 1), # Sum of absolute values. Manhattan Distance.
        l0_norm = compute_norm(x, p = 0), # Number of non-zero elements in vector
        max_channel_fraction = max(abs(x)) # Equivalent to the infinity norm
      )
    }, simplify = FALSE
  )

  # Convert Metrics To Dataframe
  df_metrics <- as.data.frame(do.call("rbind", ls_metrics))
  df_metrics[["id"]] <- colnames(mx_signatures)
  df_metrics <- df_metrics[c("id", setdiff(colnames(df_metrics), "id"))]

  # Compute versions scaled based on number of channels
  df_metrics[["shannon_index_exp_scaled"]] <- df_metrics[["shannon_index_exp"]]/n_channels
  df_metrics[["l3_norm_scaled"]] <- df_metrics[["l3_norm"]]/n_channels
  df_metrics[["l2_norm_scaled"]] <- df_metrics[["l2_norm"]]/n_channels
  df_metrics[["l1_norm_scaled"]] <- df_metrics[["l1_norm"]]/n_channels
  df_metrics[["l0_norm_scaled"]] <- df_metrics[["l0_norm"]]/n_channels

  # Add scaled
  rownames(df_metrics) <- NULL
  return(df_metrics)
}

#' Compute Pairwise Similarity or Distance Metrics Between Signatures
#'
#' Calculates a specified pairwise metric for all unique pairs of signatures in a `sigverse` collection.
#' Supported metrics are:
#' - **`cosine_similarity`**: returns values in [0, 1], where 1 indicates identical signatures
#' - **`L2`**: Euclidean distance between signatures
#' - **`L1`**: Manhattan distance between signatures
#'
#' The input may be supplied either as:
#' - A **named list** of `sigverse` signature data.frames (each validated via
#'   [sigshared::assert_signature_collection()]), or
#' - A **numeric matrix** (rows = mutation contexts, columns = signatures) with non-NULL, unique column names.
#'
#' Before computing any metric, each column is normalized to sum to 1 (via
#' `sigshared::compute_fraction()`), ensuring minor rounding errors in provided fractions do not affect results.
#'
#' @param signatures A signature collection, either as a named list of `sigverse` signature data.frames
#'   or as a numeric matrix of fractional values (contexts × signatures). Column names are used as signature IDs.
#' @param metric Character; one of `c("cosine_similarity", "L2", "L1")`.
#'   Determines which pairwise metric to compute.
#' @param format Character; one of `c("data.frame", "matrix")`.
#'   - `"data.frame"` (default): returns a long-format data.frame with columns `S1`, `S2`, and `<metric>`.
#'   - `"matrix"`: returns an N×N symmetric numeric matrix of pairwise metrics (dimnames = signature IDs).
#'     The diagonal is set to 1 for cosine similarity, and 0 for the distance metrics.
#'
#' @return Depending on `format`:
#' \describe{
#'   \item{`data.frame`}{A data.frame with one row per unique signature pair, and columns
#'     `S1`, `S2`, and the chosen metric (named after `metric`).}
#'   \item{`matrix`}{A symmetric numeric matrix of pairwise metrics, with row and column names
#'     corresponding to signature IDs.}
#' }
#'
#' @examples
#' library(sigstash)
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Long-format cosine similarities
#' sig_collection_pairwise_stats(signatures,
#'                               metric = "cosine_similarity",
#'                               format = "data.frame")
#'
#' # Symmetric matrix of Euclidean distances
#' sig_collection_pairwise_stats(signatures,
#'                               metric = "L2",
#'                               format = "matrix")
#'
#' @seealso [sig_cosine_similarity()], [sig_l2_distance()], [sig_lp_distance()], [sig_collection_stats()]
#' @export
sig_collection_pairwise_stats <- function(signatures,
                                          metric = c("cosine_similarity", "L2", "L1"),
                                          format = c("data.frame", "matrix")
                                          ){

  metric <- rlang::arg_match(metric)
  format <- rlang::arg_match(format)

  # Assert signature class and convert to matrix form, if not already
  if(is.list(signatures)){
    sigshared::assert_signature_collection(signatures)
    signatures <- sigshared::sig_collection_reformat_list_to_matrix(signatures, values = "fraction")
  }
  # If not a list, assert its a valid signature matrix (rows=channels cols=samples)
  else{
    sigshared::assert_signature_collection_matrix(signatures, must_sum_to_one = FALSE)
  }

  # Normalise any counts to fractions. If already fractions that sum to one,
  # this will make no difference.
  signatures <- apply(X = signatures, MARGIN = 2, sigshared::compute_fraction, simplify = TRUE)

  # Create pairwise combinations
  samples = colnames(signatures)
  ls_combinations = combn(x = samples, m = 2, simplify = FALSE)

  metric_function <-
    if(metric == "cosine_similarity")
      compute_cosine_similarity
  else if(metric == "L2")
      function(x, y) { compute_lp_distance(x, y, p = 2) }
  else if(metric == "L1")
    function(x, y) { compute_lp_distance(x, y, p = 1) }
  else {
    stop("Distance metric: ", metric, " has not yet been implemented. Please create issue on the sigstats github page to request this distance metric.")
  }

  # Compute pairwise metrics
  vec_metrics <- vapply(
    X = ls_combinations,
    function(combos){
      sample1 = combos[1]
      sample2 = combos[2]
      vec1 = signatures[, sample1, drop=TRUE]
      vec2 = signatures[, sample2, drop=TRUE]

      # Compute Metric
      metric_function(vec1, vec2)
    },
    FUN.VALUE = numeric(1)
  )

  # Create data.frame of sample combinations
  df_combinations <- as.data.frame(do.call("rbind", ls_combinations))
  colnames(df_combinations) <- c("S1", "S2")

  # Output matrix
  if(format == "matrix"){
    m <- matrix(NA_real_, length(samples), length(samples),
                dimnames = list(samples, samples))
    m[df_combinations$S1, df_combinations$S2] <- vec_metrics
    m[df_combinations$S2, df_combinations$S1] <- vec_metrics
    diag(m) <- if (metric == "cosine_similarity") 1 else 0
    return(m)
  }

  # Output data.frame
  df_combinations[[metric]] <- vec_metrics

  return(df_combinations)
}


# Underlying Computations -------------------------------------------------
compute_shannon_index <- function(probabilities, exponentiate = FALSE){

  # Drop zeros and avoid log(0)
  probabilities <- probabilities[probabilities > 0]

  # Compute Shannon
  shannon_index <- -sum(probabilities * log(probabilities))


  if(exponentiate){
    return(exp(shannon_index))
  }

  return(shannon_index)
}

compute_kl_divergence <- function(p, q, pseudocount = 1e-12, base = exp(1)){
  if(all(p == 0)) return(0)
  p <- p / sum(p)
  sum(p*log(p/(q+pseudocount), base = base), na.rm = TRUE)
}

compute_gini <- function(x) {
  x <- sort(x)
  n <- length(x)
  if (all(x == 0)) return(0)  # define Gini = 0 when all entries are 0
  G <- sum((2 * seq_len(n) - n - 1) * x)
  G <- G / (n^2 * mean(x))
  G * n/(n-1)

}

compute_norm <- function(vec, p){
  if(is.infinite(p) & p>0){
    return(max(abs(vec)))
  }

  if(p==0){
    return(sum(vec != 0))
  }

  (sum(abs(vec)^p))^(1/p)
}

# Scale: (boolean) scale to the length of vec1 & vec2
compute_lp_distance <- function(vec1, vec2, p){
  vec <- vec1 - vec2
  compute_norm(vec, p = p)
}




compute_cosine_similarity <- function(x, y){
  as.numeric(lsa::cosine(x, y))
}




# Bootstrap Stats ---------------------------------------------------------

#' Calculate Experimental P-Value Based on Bootstrap Contributions
#'
#' For a single signature, calculate the proportion of bootstraps that have a percentage contribution below a specified threshold. This value represents an 'experimental' P-value. For example, if P < 0.05, there is a <5% chance of the signature having a contribution >= \code{threshold} by chance.
#'
#' @param bootstraps A numeric vector describing the percentage contributions of a single signature in different bootstraps.
#' @param threshold A numeric value specifying the minimum contribution threshold (typically 0.05). This is referred to as the sparsity threshold.
#'
#' @return The proportion of bootstraps with contributions below \code{threshold}. This represents an experimental p-value (referred to as the sparsity p-value).
#' @export
#'
#' @examples
#' bootstraps <- c(0.01, 0.03, 0.07, 0.02, 0.05)
#' threshold <- 0.05
#' sig_compute_experimental_p_value(bootstraps, threshold)
sig_compute_experimental_p_value <- function(bootstraps, threshold) {
  sum(bootstraps < threshold) / length(bootstraps)
}


#' Summarise Signature Bootstraps
#'
#' This function computes summary statistics for signature-level bootstraps,
#' including quantiles, min/max contribution proportions, and an experimental p-value.
#'
#' @param bootstraps A data frame containing bootstrap results in a sigverse-style format.
#'                   See [sigshared::example_bootstraps()] for details.
#' @param threshold The minimum proportion of mutations that must be explained by a signature for it to be
#'                  considered present in a bootstrap. This threshold is used for calculating
#'                  the experimental p-value via [sigstats::sig_compute_experimental_p_value()].
#' @inheritParams boxplotstats::calculate_boxplot_stats
#'
#' @return A data frame containing the following columns:
#' \item{signatures}{Names of the signatures.}
#' \item{quantiles, iqr, min, max,outlier_thresholds}{Statistical summary columns (quantiles, min, and max) of the contribution proportions across bootstraps.}
#' \item{outliers}{either a list column or a character of '|' outlier strings (depending on \code{outliers_as_strings} argument}
#' \item{p_value}{Computed experimental p-value for each signature, based on the threshold.}
#' @export
#'
#' @examples
#' library(sigshared)
#' bootstraps <- example_bootstraps()
#' sig_summarise_bootstraps(bootstraps, threshold = 0.05)
sig_summarise_bootstraps <- function(bootstraps, threshold = 0.05, outliers_as_strings = FALSE) {

  # Ensure that the input data frame is valid for bootstraps
  sigshared::assert_bootstraps(bootstraps)

  # Apply a function over each signature group, computing summary stats and p-value
  bootstrap_summary <- tapply(X = bootstraps$contribution,
                              INDEX = bootstraps$signature,
                              FUN = function(contributions) {

                                # Calculate statistical summary (boxplot stats) for the contributions
                                summary_stats <- boxplotstats::calculate_boxplot_stats(contributions, return_dataframe = TRUE, outliers_as_strings = outliers_as_strings)

                                # Compute the experimental p-value for each signature based on the threshold
                                summary_stats$p_value <- sigstats::sig_compute_experimental_p_value(contributions, threshold)

                                return(summary_stats)
                              }, simplify = FALSE)

  # Combine the list of summary data frames into a single data frame
  result_df <- do.call(rbind, bootstrap_summary)

  # Move signature names from row names to a new 'signatures' column
  result_df$signatures <- rownames(result_df)
  rownames(result_df) <- NULL

  # Reorder columns to have 'signatures' as the first column
  result_df <- result_df[c("signatures", setdiff(colnames(result_df), "signatures"))]

  return(result_df)
}
