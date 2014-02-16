\name{vec.to.matrix}
\alias{vec.to.matrix}
\title{Convert a by-rows vector to a genotype matrix}
\usage{
vec.to.matrix(cv)
}
\arguments{
  \item{cv}{vector containing \code{k(k+1)/2} genotype
  counts. All non-negative integers.}
}
\value{
a \code{k x k} symmetrical matrix of genotype counts.
}
\description{
If there are \code{k} alleles, there are \code{k(k+1)/2}
possible genotypes. Convert a row- wise vector of genotype
counts into a \code{k x k} symmetrical matrix of genotype
counts. Genotype counts should be in the order: \code{a11,
a21, a22, a31, a32, ..., a}
}
\examples{
vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
}

