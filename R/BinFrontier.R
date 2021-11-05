BinFrontier <- function(treatment, dataset, formula, metric, QOI, breaks, match.on, verbose) {

  treat <- dataset[[treatment]]

  if (verbose) cat("Computing bins...\n")

  if (endsWith(metric, "median")) {
    bins.list <- getBinsAtMedian(dataset, match.on, treat, metric)
  }
  else {
    bins.list <- getBins(dataset, match.on, breaks)
  }
  strataholder <- assignToBins(dataset, match.on, bins.list)

  if (verbose) cat("Calculating frontier...\n")

  if (QOI == "FSATE") {
    frontier <- binsToFrontierFSATE(strataholder, treat, metric)
  }
  else if (QOI == "SATT") {
    frontier <- binsToFrontierSATT(strataholder, treat, metric)
  }

  if (verbose) cat("Done!\n")

  return(frontier)
}
