\name{genotypeList.to.matrix}
\alias{genotypeList.to.matrix}
\title{Convert a list of genotypes into a genotype count matrix}
\usage{
genotypeList.to.matrix(g)
}
\arguments{
  \item{g}{List of text objects indicating genotypes.
  Alleles are separated by \dQuote{/}}
}
\value{
matrix of \eqn{k x k} genotype counts
}
\description{
Genotype lists as are used by packages `genetics` and
`adegenet` are converted to an array of genotype counts.
This function requires package `genetics`
}
\examples{
g <- c(rep("a/b",3),
		"b/b",
		rep("a/c", 5),
		rep("b/c", 18),
		"c/c",
		rep("a/d",3),
		rep("b/d", 7),
		rep("c/d", 5),
		rep("d/d", 2))
genotypeList.to.matrix(g)
}

