#' Add Signatures to a Combined Signature Model
#'
#' This function takes a signature collection and a model as input and adds selected signatures to a combined signature model. The output is a data.frame in the 'combined_signature_model' style data.frame format from the Sigverse package.
#'
#' @param signatures A signature collection, typically a list of data.frames where each data.frame represents a signature with specific attributes.
#' @param model A numeric vector representing the contribution of each signature to the combined model. The sum of the values in this vector should be less than or equal to 1.
#'
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
sig_combine <- function(signatures, model){

  # Assertions
  sigshared::assert_signature_collection(signatures)
  assertions::assert_numeric_vector(model)

  model_signatures <- names(model)
  assertions::assert_subset(model_signatures, names(signatures))
  assertions::assert(sum(model) <= 1, msg = 'Contributions of all signatures in model should add up to <= 1, not [{sum(model)}]')

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

  # Perform the addition
  return(df_signatures_combined)
}
