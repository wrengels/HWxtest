\name{xcountCutoff}
\alias{xcountCutoff}
\title{Determine immediately whether number of tables is over a limit}
\usage{
\method{generic}{xcountCutoff}(x, ...)
}
\arguments{
  \item{m}{vector containing the numbers of alleles of each
  type. Length must be at least 2 and all must be positive
  integers. It can also be a matrix of genotype counts, but
  not a vector of genotype counts.}

  \item{...}{parameters passed to \code{xcountCutoff}}

  \item{cutoff}{Is the number of tables above or below this
  value?}
}
\value{
TRUE or FALSE depending on whether the table count is above
or below \code{cutoff}
}
\description{
Calling \code{scountCutoff} gives you a quick answer to
whether the number of tables is over a given cutoff. It is
useful in deciding whether to analyze a data set with
\code{\link{xtest}} or \code{\link{mtest}}. This function
is used by \code{\link{hw.test}} and not normally called
directly by the user.
}
\examples{
#
alleles <- c(15, 14, 11, 12, 2, 2, 1, 3)
if(xcountCutoff(alleles)) cat("There are too many tables")
}

