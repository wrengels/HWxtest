\name{show.hwtest}
\alias{show.hwtest}
\title{Display results from \code{\link{hw.test}}}
\usage{
show.hwtest(h, detail = NA, statName = NA)
}
\arguments{
  \item{h}{the results of a call to \code{\link{hw.test}}}

  \item{detail}{may be used to change the detail of the
  output}

  \item{statName}{may be used to specify which statistic to
  use if \code{detail} is 1.}
}
\value{
none
}
\description{
\code{show.hwtest} is called by \code{print}, but it can
also be called directly by the user if the object is a list
of \code{hwtest} objects, or even a list of such lists.
}

