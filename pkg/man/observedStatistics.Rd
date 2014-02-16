\name{observedProb}
\alias{observedLLR}
\alias{observedProb}
\alias{observedU}
\alias{observedX2}
\title{Compute observed statistics for a genotype count matrix}
\usage{
observedProb(c)

observedLLR(c)

observedU(c)

observedX2(c)
}
\arguments{
  \item{c}{Matrix of observed genotype counts. Each number
  should be a non-negative integer, and matrix is \eqn{k x
  k}.}
}
\value{
the observed statistic
}
\description{
Four measures of fit to Hardy-Weinberg for a given set of
genotype counts may be computed.
}
\details{
* \code{observedProb} The probability of the observed set
under the HW null and with the allele counts fixed.

* \code{observedLLR} The log-likelihood ratio of the
observed set

* \code{observedU} The observed U-score. Positive values
indicate an excess of homozygotes and negative ones imply
too many heterozygotes

* \code{observedX2} The classical \dQuote{chi-squared}
statistic
}
\examples{
t <- vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
observedStats <- c(observedProb(t), observedLLR(t), observedU(t), observedX2(t))
}

