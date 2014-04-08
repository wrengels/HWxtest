\name{genepop.to.genind}
\alias{genepop.to.genind}
\title{Imports a \code{.txt} file in \code{GenePop} format into an object of type \code{genind}}
\usage{
genepop.to.genind(name)
}
\arguments{
  \item{name}{the name of a file in \code{GenePop} format}
}
\value{
an object of class \code{genind}
}
\description{
The main work is done by the function
\code{adegenet::import2genind}. However, that function
requires text files with an extension of \code{.gen},
whereas such files usually have extension \code{.txt}. The
sole purpose of this function is to work around the
\dQuote{.gen} requirement.
}

