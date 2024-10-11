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
