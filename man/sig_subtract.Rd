% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_maths.R
\name{sig_subtract}
\alias{sig_subtract}
\alias{\%-\%}
\title{Subtract Signatures/Catalogues}
\usage{
sig_subtract(signature1, signature2)

signature1 \%-\% signature2
}
\arguments{
\item{signature1, signature2}{sigverse signature/catalogue data.frames. See \code{\link[sigshared:signature]{sigshared::example_signature()}} or \code{\link[sigshared:catalogue]{sigshared::example_catalogue()}}.}
}
\value{
a data.frame representing a sigverse signature.
}
\description{
Subtracts signature2 from signature1 and returns result.
}
\examples{
library(sigstash)

# Load a signature collection
signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")

# Subtract signatures
signatures[['SBS3']] \%-\% signatures[['SBS4']]

# Identical approach using full function name
sig_subtract(signatures[['SBS3']], signatures[['SBS4']])

}
