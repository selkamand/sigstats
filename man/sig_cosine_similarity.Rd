% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_maths.R
\name{sig_cosine_similarity}
\alias{sig_cosine_similarity}
\title{Calculate Cosine Matrix Between Two Signatures}
\usage{
sig_cosine_similarity(signature1, signature2, assume_sensible_input = FALSE)
}
\arguments{
\item{signature1, signature2}{sigverse signature data.frames}

\item{assume_sensible_input}{to drastically speed up the similarity function,
we an simply assume that both inputs are valid signature objects, and the channels are sorted.
If speed is essential, perform these checks upstream (flag)}
}
\value{
a number between 0 and 1 representing cosine similarity
}
\description{
Computes cosine similarity between each pair of signatures in a sigverse signature collection
}
\examples{
library(sigstash)

# Load a signature collection
signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")

# Compute cosine similarity between two signatures
sig_cosine_similarity(signatures[["SBS1"]], signatures[["SBS2"]])

}
