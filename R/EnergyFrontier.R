EnergyFrontier <- function(treatment, dataset, formula, metric, QOI, distance.mat = NULL, verbose){

  treat <- dataset[[treatment]]

  if (verbose && is.null(distance.mat)) cat("Computing distance matrix...\n")

  if (is.null(distance.mat)) {
    covs.mat <- scale(get.covs.matrix(formula, dataset))
    distance.mat  <- as.matrix(dist(covs.mat))
  }
  else {
    if (!is.matrix(distance.mat) || !is.numeric(distance.mat) ||
        nrow(distance.mat) != length(treat) || ncol(distance.mat) != length(treat) ||
        anyNA(distance.mat)) {
      customStop("'distance.mat' must be NULL or an N x N numeric matrix when metric = \"Energy\".", "makeFrontier()")
    }
  }

  if (verbose) cat("Calculating frontier...\n")

  if (QOI == "SATE") {
    frontier <- energyToFrontierSATE(distance.mat, treat)
  }
  else if (QOI == "FSATE") {
    frontier <- energyToFrontierFSATE(distance.mat, treat)
  }
  else if (QOI == "SATT") {
    frontier <- energyToFrontierSATT(distance.mat, treat)
  }

  if (verbose) cat("Done!\n")

  return(frontier)
}
