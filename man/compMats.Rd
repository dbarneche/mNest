\name{compMats}
\alias{compMats}
\title{Performs a comparison between matrices via a given function}
\usage{
  compMats(x, func, mat1, mat2)
}
\arguments{
  \item{x}{index of mat2}

  \item{mat1}{original matrix}

  \item{mat2}{matrix to be compared}

  \item{func}{function used to compare}
}
\value{
  lgical. If TRUE then matrices are identical according to
  the parameter requested
}
\description{
  Are all outputs from a given function the same between
  two matrices?
}

