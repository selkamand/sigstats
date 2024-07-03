`%||%` <- function(lhs, rhs){
 if(is.null(lhs))
   return(rhs)
  else
    return(lhs)
}

#' Add Signatures to a Combined Signature Model
#'
#' This function takes a signature collection and a model as input and adds selected signatures to a combined signature model. The default output is a data.frame in the 'combined_signature_model' style data.frame format from the Sigverse package.
#'
#' @param signatures A signature collection, typically a list of data.frames where each data.frame represents a signature with specific attributes.
#' @param model A named numeric vector representing the contribution of each signature to the combined model.
#' The names correspond to the signatures, and the values represent their contributions.
#' The sum of the values in this vector should be less than or equal to 1.
#' @param format A character string indicating the output format.
#' If "combined", the function returns a 'combined_signature_model' data.frame where each row represents a contribution for a particular channel from a single signature (duplicate channels are not collapsed).
#' If "signature", the function returns the data in the sigverse signature format, representing a novel signature created by combining the signatures in the collection according to the ratios described by the model.
#' @param verbose enables detailed output messages (flag).
#' @return A data.frame in the 'combined_signature_model' style containing the selected signatures and their modified fractions based on the model.
#'
#' @seealso \code{\link{assert_signature_collection}}, \code{\link{assert_numeric_vector}}, \code{\link{assert_subset}}
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
#' # Visualise using sigvis
#'
#' @export
sig_combine <- function(signatures, model, format = c("signature", "combined"), verbose=FALSE){

  # Settings
  tolerance <- 5e-07 # Allows models to sum to anything < 1.000001

  # Null replacements
  model <- model %||% numeric(0)

  # Assertions
  sigshared::assert_signature_collection(signatures)
  assertions::assert_numeric_vector(model)
  format <- rlang::arg_match(format)

  model_signatures <- names(model)
  assertions::assert_subset(model_signatures, names(signatures))
  assertions::assert(sum(model) <= 1 + tolerance, msg = 'Contributions of all signatures in model should add up to <= 1, not [{sum(model)}]')
  assertions::assert_no_duplicates(model_signatures)

  # Deal with the case that model vector has no length.
  # If format="signature" we can just return an empty signature (all fraction values = 0)
  # If the format="combined" (individual signature contributions preserved in dataframe) there is no
  # Sensible way to represent an model with no contributions, and so an error is thrown
  if(is.null(model) | length(model_signatures) == 0){
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
  ls_signatures <- lapply(model_signatures, \(signame){
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
  stats::aggregate(data = signature_combination, fraction ~ type + channel, FUN = \(frac){ sum(frac)})
}


#' Calculate Cosine Matrix Between Two Signatures
#'
#' Computes cosine similarity between each pair of signatures in a sigverse signature collection
#'
#' @param signature1,signature2 sigverse signature data.frames
#' @param assume_sensible_input to drastically speed up the similarity function,
#' we an simply assume that both inputs are valid signature objects, and the channels are sorted.
#' If speed is essential, perform these checks upstream (flag)
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
    sig1_type_channel_id = paste(signature1[['type']], signature1[['channel']])
    sig2_type_channel_id = paste(signature2[['type']], signature2[['channel']])
    assertions::assert_equal(sig1_type_channel_id, sig2_type_channel_id, msg = 'Can NOT calculate cosine similarity for two signatures/catalogues which have different types or channels.')
  }
  cosine <- sim_cosine(signature1[['fraction']], signature2[['fraction']])

  # If NaN (e.g.if all fractions = 0) because replace with 0
  if(is.nan(cosine)) cosine <- 0

  return(cosine)
}


sim_cosine <- function(x, y){
  as.numeric(lsa::cosine(x, y))
}

#' Reconstruct a catalogue from a signature
#'
#' Generates reconstructed substitution profile (catalogues) generated by proportionally combining each estimated contribution.
#'
#' @param signature a data.frame in sigverse signature format
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


#' Subtract Signature
#'
#' Subtracts signature2 from signature1 (factions) and returns result.
#' Works for only signatures. To translate to a catalogue see [sig_reconstruct()]
#'
#' @param signature1,signature2 sigverse signature data.frames
#'
#' @return a data.frame representing a sigverse signature
#' @export
#'
#' @examples
#' library(sigstash)
#'
#' # Load a signature collection
#' signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Subtract signatures
#' sig_subtract(signatures[['SBS3']], signatures[['SBS4']])
#'
sig_subtract <- function(signature1, signature2){
  sigshared::assert_signature(signature1, must_sum_to_one = FALSE)
  sigshared::assert_signature(signature2, must_sum_to_one = FALSE)

  # Ensure signatures are sorted the same way 'type channel' match 1:1 (including order)
  sig1_type_channel_id = paste(signature1[['type']], signature1[['channel']])
  sig1_type_channel_id = paste(signature2[['type']], signature2[['channel']])
  assertions::assert_equal(sig1_type_channel_id, sig1_type_channel_id)

  sig_result <- signature1
  sig_result[['fraction']] <- signature1[['fraction']] - signature2[['fraction']]

  return(sig_result)
}
