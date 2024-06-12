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
