\name{nullBinMats}
\alias{nullBinMats}
\title{Generates null binary matrices following the null model 2 in Bascompte et al. 2003 PNAS.
Written by M.M.Pires, optmized by Barneche (in progress)}
\usage{
  nullBinMats(mat, iter = 100)
}
\arguments{
  \item{mat}{Binary matrix}

  \item{iter}{Number of iterations (recommended 1000)}
}
\value{
  mat.t array of randomized matrices
}
\description{
  This null created random matrices by randomly resorting
  the 1's among the matrix cells according to marginal
  totals of rows and columns.
}
\details{
  Each cell has a probability of being filled that is
  proportional to the number of occurrences of individuals
  in sites: cij = 1/2*(Pi/C + Pj/R) where Pi= row sums; Pj
  = column sums; C = number of columns; and R = number of
  rows.
}
\references{
  "Bascompte J, Jordano P, Melian CJ and Olesen JM (2003)
  The nested assembly of plant-animal mutualistic networks.
  PNAS, 100(16):9383-9387"
}

