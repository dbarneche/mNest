\name{nullMatsVegan}
\alias{nullMatsVegan}
\title{Generate random null matrices}
\usage{
  nullMatsVegan(m, iter, model, ...)
}
\arguments{
  \item{m}{A quantitative matrix}

  \item{iter}{Number of iteractions, i.e. number of random
  matrices to be created}

  \item{model}{Function to be chosen. See \code{Details}.}

  \item{...}{Further arguments from \code{permatswap} or
  \code{permatfull}}
}
\value{
  a list with iter random matrices
}
\description{
  Uses algorithms from package \code{vegan} to generate
  random matrices
}
\details{
  If \code{model} = 1, algorithm from \code{permatswap} is
  called. If \code{model} = 2, algorithm from
  \code{permatfull} is called.
}

