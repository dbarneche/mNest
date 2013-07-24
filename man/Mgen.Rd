\name{Mgen}
\alias{Mgen}
\title{Mgen algorithm (modified by Barneche on 2013)}
\usage{
  Mgen(m, n = sum(m), zs = FALSE, rep.cell = TRUE)
}
\arguments{
  \item{m}{quantitative matrix}

  \item{n}{total number of occurrences in the matrix (if
  individuals, default is the total sum of the matrix)}

  \item{zs}{logical. See \code{Details}}

  \item{rep.cell}{logical. See \code{Details}}
}
\value{
  a randomized matrix
}
\description{
  This is a generic function to build null models.
}
\details{
  \code{Mgen} is general in the sense that it allows any
  type of probability matrix to be used for constructing
  the simulated matrices. It does not, however, constrain
  rows and column totals, nor does it constrain
  connectance. If rep.cell=TRUE, repeated interactions are
  added, thus generating a quantitative matrix with cell
  values as positive integers.  If rep.cell=FALSE, no
  repeated assignment of interactions is allowed, thus
  generating a binary matrix of zeros and ones. Please note
  that when rep.cell=FALSE the number of interactions to be
  assinged must be equal or lower than the number of cells
  in the matrix. if \code{zs} is FALSE it does not allow
  any column or row to be empty, i.e. sum up zero
  (default); if TRUE, it does. If \code{rep.cell} is TRUE
  it returns a quantitative matrix (default); if FALSE, it
  returns a binary matrix.
}
\references{
  "Vázquez DP, Chacoff N and Cagnolo L (2009) Evaluating
  multiple determinants of the structure of mutualistic
  networks. Ecology, 90:2039-2046"
}

