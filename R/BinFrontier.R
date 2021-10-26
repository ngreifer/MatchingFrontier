BinFrontier <- function(treatment, dataset, formula, metric, QOI, breaks, match.on, verbose) {

  treat <- dataset[[treatment]]

  if (verbose) cat("Computing bins...\n")

  if (endsWith(metric, "median")) {
    binnings <- getBinsAtMedian(dataset, match.on, treat, metric)
  }
  else {
    binnings <- getBins(dataset, match.on, breaks)
  }

  if (verbose) cat("Calculating frontier...\n")

  if (QOI == "FSATE") {
    frontier <- binsToFrontierFSATE(binnings, treat, metric)
  }
  else if (QOI == "SATT") {
    frontier <- binsToFrontierSATT(binnings, treat, metric)
  }

  if (verbose) cat("Done!\n")

  return(frontier)
}
