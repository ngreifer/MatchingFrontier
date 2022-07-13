BinFrontier <- function(treatment, dataset, formula, metric, QOI, match.on, keep.n.equal, verbose) {

  treat <- dataset[[treatment]]

  if (is.null(attr(metric, "breaks"))) breaks <- "sturges"
  else breaks <- attr(metric, "breaks")

  if (verbose) cat("Computing bins...\n")

  if (is.character(breaks) && identical(tolower(breaks), "median")) {
    bins.list <- getBinsAtMedian(dataset, match.on, treat, metric)
  }
  else {
    bins.list <- getBins(dataset, match.on, breaks)
  }
  strata <- assignToBins(dataset, match.on, bins.list)

  if (verbose) cat("Calculating frontier...\n")

  if (QOI == "FSATE") {
    frontier <- binsToFrontierFSATE(strata, treat, metric, verbose, keep.n.equal)
  }
  else if (QOI == "SATT") {
    frontier <- binsToFrontierSATT(strata, treat, metric, verbose, keep.n.equal)
  }

  if (verbose) cat("Done!\n")

  return(frontier)
}
