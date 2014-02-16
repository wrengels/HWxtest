\name{fillUpper}
\alias{fillUpper}
\title{Make a matrix symmetrical, using only the lower-left half}
\usage{
fillUpper(t)
}
\arguments{
  \item{t}{a matrix of non-negative integers representing
  genotype counts. At least 2x2}
}
\value{
a matrix of the same dimensions as \code{t} but
symmetrical. Upper-right is replaced.
}
\description{
In a matrix of genotype counts, \code{a[i,j]} and
\code{a[j,i]} both represent the same heterozygote if
\code{i!=j}. So only half the matrix need be used. If
\code{i==j} then it is a homozygote.
}
\examples{
d <- matrix(1:25,5,5)
fillUpper(d)
}

