\name{matrix.to.vec}
\alias{matrix.to.vec}
\title{Extract a by-rows vector from a genotype matrix}
\usage{
matrix.to.vec(t)
}
\arguments{
  \item{t}{square matrix of genotype counts. Only the
  lower-left half is used.}
}
\value{
a vector of length \eqn{k(k-1)/2} containing the genotype
counts
}
\description{
Given a square matrix of size \eqn{k x k}, generate a
vector of genotype counts of size \eqn{k(k-1)/2}
}
\examples{
t <- vec.to.matrix(c(0,3,1,5,18,1,3,7,5,2))
v <- matrix.to.vec(t)
}

