\name{plotRarefac}
\alias{plotRarefac}
\title{Plots rarefaction curves with 95 percent CI}
\usage{
  plotRarefac(data, x, ID, permut)
}
\arguments{
  \item{data}{Central database}

  \item{x}{an identity element}

  \item{ID}{identity vector to be matched with x}

  \item{permut}{Number of iterations}
}
\value{
  individual plots for each site via \code{quartz}
}
\description{
  Calculates sample-based rarefaction curves via function
  \code{specaccum} from package \code{vegan}.
}

