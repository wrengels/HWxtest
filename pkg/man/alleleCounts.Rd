\name{alleleCounts}
\alias{alleleCounts}
\title{Extract allele counts from matrix of genotypes}
\usage{
alleleCounts(g)
}
\arguments{
  \item{g}{a matrix of non-negative integers representing
  genotype counts}
}
\value{
a vector of allele counts
}
\description{
Create a vector of the number of each allele. The sum will
be double that of genotpes.
}
\examples{
d <- matrix(1:25,5,5)
alleleCounts(d)
}

