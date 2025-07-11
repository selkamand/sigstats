`%||%` <- function(lhs, rhs){
 if(is.null(lhs))
   return(rhs)
  else
    return(lhs)
}


# Operations --------------------------------------------------------------


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

