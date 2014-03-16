\name{print.hwtest}
\alias{print.hwtest}
\title{S3 Method for printing \code{hwtest} objects}
\usage{
\method{print}{hwtest}(h, detail = h$detail, ...)
}
\arguments{
  \item{h}{the results from a call to
  \code{\link{hw.test}}}

  \item{detail}{If unspecified, the level of detail used in
  the \code{hw.test} call will be applied.}

  \item{...}{other parameters passed to \code{print}}
}
\description{
Prints test results (\code{hwtest}) objects depending on
how much detail is provided.
}

