\name{reportNestedness}
\alias{reportNestedness}
\title{Reports observed WNODF and NODF against 95%CI from null models.}
\usage{
  reportNestedness(bin, nullbin, quant, nullquant, output=TRUE, scales, addInfo)
}
\arguments{
  \item{bin}{A list of observed binary NODF summaries}
  \item{nullbin}{A list of random-generated binary NODF summaries}
  \item{quant}{A list of observed quantitative WNODF summaries}
  \item{nullquant}{A list of observed quantitative WNODF summaries}
  \item{output}{should results be written to output folder? (default=TRUE)}
  \item{scales}{if \code{output} is TRUE, scales is incorporated into ouput name. 
                Must be either 'site' or 'samples'
                }
  \item{addInfo}{If \code{output} is TRUE and \code{scales} is 'samples',
                 addInfo is incorporated into ouput name.
                 }
}
\value{
  a table of summaries comparing observed mean against 95%CI
}
\description{
  Reports observed WNODF and NODF against 95%CI from null models.
}
