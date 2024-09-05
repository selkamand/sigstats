#' Convert Signature Model to String Representation
#'
#' This function converts a named numeric vector model representing signature contributions to a formatted string.
#' By default, contributions are displayed as percentages, with options to customize formatting.
#' Unexplained portions of the model (i.e., contributions not accounted for by the signatures) can be optionally displayed.
#'
#' @inheritParams sig_combine
#' @param signature_sep A character string used to separate each signature name from its value (default is " = ").
#' @param pair_sep A character string used to separate different signatures (default is " | ").
#' @param show_percent A logical indicating whether to format values as percentages (default is TRUE).
#' @param significant_digits The number of significant digits to display for the contributions (default is 2).
#' @param include_unexplained A logical flag indicating whether to include the unexplained portion of the model (default is TRUE).
#' @param unexplained_label A character string to label the unexplained portion of the model (default is "?").
#' @return A string representing the signature model with contributions formatted according to the specified options.
#'
#' @examples
#' model <- c(SBS1 = 0.4, SBS2 = 0.6)
#' sig_model_to_string(model)
#'
#' # Example with no unexplained portion, significant digits and percentages
#' sig_model_to_string(
#'   model,
#'   significant_digits = 3,
#'   show_percent = FALSE,
#'   include_unexplained = FALSE
#'   )
#'
#' @export
sig_model_to_string <- function(model,
                                signature_sep = " = ",
                                pair_sep = " | ",
                                show_percent = TRUE,
                                significant_digits = 2,
                                include_unexplained = TRUE,
                                unexplained_label = "?") {

  # Ensure model is non-null, defaulting to an empty numeric vector
  model <- model %||% numeric(0)

  # Return an empty string if the model is empty
  if (length(model) == 0) {
    return("")
  }

  # Validate that the input is a named numeric vector
  assertions::assert_numeric_vector(model)

  # Ensure the model has names corresponding to the signatures
  if (length(model) > 0) {
    assertions::assert(
      !is.null(names(model)),
      msg = "model must be a named numeric vector where names correspond to signatures, and values represent the proportion of mutations the signature explains. No names could be found."
    )
  }

  # Sort the model in descending order by contribution value
  model <- sort(model, decreasing = TRUE)

  # Add unexplained portion to the model if requested
  if (include_unexplained) {
    unexplained_contribution <- 1 - sum(model)
    unexplained <- c(unexplained_contribution)
    names(unexplained) <- unexplained_label
    model <- c(model, unexplained)
  }

  # Format values either as percentages or significant digits based on the 'show_percent' flag
  if (show_percent) {
    model <- fmt_percent(model, digits = significant_digits)
  } else {
    model <- signif(model, digits = significant_digits)
  }

  # Combine the formatted model into a string with the specified separators
  paste0(paste0(names(model), signature_sep, model), collapse = pair_sep)
}

# Helper function to format numeric values as percentages
# This function takes a numeric vector, multiplies by 100, and formats the result with the specified number of significant digits.
# @param values A numeric vector to be converted to percentages.
# @param digits The number of significant digits to include in the formatted output (default is 3).
fmt_percent <- function(values, digits = 3) {
  # Retain original names of the input vector
  original_names <- names(values)

  # Convert numeric values to percentage strings with the specified number of digits
  formatted_values <- paste0(signif(100 * values, digits = digits), "%")

  # Restore the original names to the formatted vector
  names(formatted_values) <- original_names

  return(formatted_values)
}
