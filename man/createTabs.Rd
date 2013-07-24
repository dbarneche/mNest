\name{createTabs}
\alias{createTabs}
\title{Creates tables of traits (rows) and scale (columns)}
\usage{
  createTabs(data, vals, target, scale, func,
    assign = FALSE, tabname)
}
\arguments{
  \item{data}{Central database}

  \item{vals}{Values to be used in table}

  \item{target}{Traits to be used}

  \item{scale}{Scale to be used for columns}

  \item{func}{Function to pe parsed to \code{tapply}}

  \item{assign}{Should final table be assigned to an
  object?}

  \item{tabname}{Name to be assigned to final table}
}
\value{
  assigns tabname a table of target (rows) x scale
  (columns) taking vals as input
}
\description{
  Creates tables of traits (rows) and scale (columns)
}

