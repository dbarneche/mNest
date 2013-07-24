\name{getTraits}
\alias{getTraits}
\title{Adds traits to the central dataset}
\usage{
  getTraits(data, traits, matcH, targets, verbose = FALSE)
}
\arguments{
  \item{data}{Central database}

  \item{traits}{Traits database}

  \item{matcH}{Name of column or vector from which to
  conduct the merge}

  \item{targets}{Vector of columns names to be merged}

  \item{verbose}{If TRUE, print messages to screen, good
  for isolating problems}
}
\value{
  Central database with additional traits columns
}
\description{
  Adds traits to the central dataset
}

