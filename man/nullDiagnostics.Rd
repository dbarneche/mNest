\name{nullDiagnostics}
\alias{nullDiagnostics}
\title{Performs dignostics between original and null matrices}
\usage{
  nullDiagnostics(original, output, iters, func,
    verbose = TRUE)
}
\arguments{
  \item{original}{Original matrix}

  \item{output}{List containing all random matrices
  generated via \code{nullMatsVegan}}

  \item{iters}{Index of matrices to be compared (default
  uses the full length of output, i.e. original number of
  iterations)}

  \item{func}{Function used to compare matrices via
  \code{compMats}}

  \item{verbose}{If TRUE, print messages to screen, good
  for isolating problems}
}
\value{
  Output messages from \code{compMats} via \code{sapply}
}
\description{
  Performs dignostics for null models 1-4
}

