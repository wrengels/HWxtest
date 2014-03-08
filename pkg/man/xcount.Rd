\name{xcount}
\alias{xcount}
\title{Find Exact Number of Genotype Tables}
\usage{
xcount(x, ...)
}
\arguments{
  \item{m}{vector containing the numbers of alleles of each
  type. Length must be at least 2 and all must be positive
  integers.}

  \item{safety}{Stop execution if the approximate table
  number obtained from \code{acount()} is more than this
  cutoff.}

  \item{safeSecs}{Time limit in seconds. Another safety
  feature to prevent getting stuck in a too-long
  computation}
}
\value{
The exact number of tables
}
\description{
Use \code{xcount} to determine the exact number of tables
(i.e., genotype numbers) for a given set of allele counts.
This method enumerates all tables, and is best when the
total number is less than 10^10 or so.
}
\examples{
# Allele counts from human Rh locus. Guo and Thompson, 1992, Figure 1
#
alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
xcount(alleles)
}

