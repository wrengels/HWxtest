\name{acount}
\alias{acount}
\title{Find Approximate Number of Genotype Tables}
\usage{
acount(m, ...)
}
\arguments{
  \item{m}{vector containing the numbers of alleles of each
  type. Length must be at least 2. All items are positive
  integers.}
}
\value{
The approximate number of tables.
}
\description{
Use \code{acount} to obtain the approximate number of
genotype tables for a given set of allele counts. This
method uses a normal approximation and is much faster than
enumerating the tables.
}
\examples{
# Allele counts from human Rh locus. Guo and Thompson, 1992, Figure 1
#
alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
acount(alleles)
# This approximation may be compared with the exact value of 250552020
#
ld <- c(6329, 319, 47, 2773, 75, 6702, 14, 2, 333)
acount(ld)
#
# This is an example where the number of tables is too large for a full enumeration.
}

