source("R/import.R")

dataDir  <-  "datasets"
FigsDir  <-  "figures"
OutsDir  <-  "output"

packages  <-  c("vegan", "bipartite")
packages  <-  invisible(sapply(packages, isPack, verbose=TRUE))
rm(packages)
