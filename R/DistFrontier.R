# DistFrontier <- function(treatment, dataset, metric, QOI, match.on, distance.mat = NULL){
DistFrontier <- function(treatment, dataset, formula, metric, QOI, distance.mat = NULL, verbose){

  treat <- dataset[[treatment]]

  if (verbose && metric != "Custom") cat("Computing distance matrix...\n")

  if (metric == "Mahal") {
    covs.mat <- get.covs.matrix(formula, dataset)
    distance.mat <- calculateMdist(covs.mat, treat)
  }
  else if (metric == "Euclid") {
    covs.mat <- get.covs.matrix(formula, dataset)
    distance.mat <- calculateEdist(covs.mat, treat)
  }
  else if (metric == "Custom") {
    if (is.null(distance.mat)) {
      customStop("'distance.mat' must be specified when metric = \"Custom\".", "makeFrontier()")
    }
    if (!is.matrix(distance.mat) || !is.numeric(distance.mat) ||
        nrow(distance.mat) != sum(treat == 1) || ncol(distance.mat) != sum(treat == 0) ||
        anyNA(distance.mat)) {
      customStop("'distance.mat' must be an N1 x N0 numeric matrix when metric = \"Custom\".", "makeFrontier()")
    }
  }

  if (verbose) cat("Calculating frontier...\n")

  if (QOI == "FSATE") {
    frontier <- distToFrontierFSATE(distance.mat, treat)
  }
  else if (QOI == "FSATT") {
    frontier <- distToFrontierFSATT(distance.mat, treat)
  }

  if (verbose) cat("Done!\n")

  return(frontier)
}
