DistFrontier <- function(treatment, dataset, formula, metric, QOI, verbose){

  treat <- dataset[[treatment]]

  if (is.null(attr(metric, "distance.mat"))) distance.mat <- "mahalanobis"
  else distance.mat <- attr(metric, "distance.mat")

  if (is.character(distance.mat) && length(distance.mat) == 1L && distance.mat != "custom") {
    if (verbose) cat("Computing distance matrix...\n")

    distance.mat <- match_arg(distance.mat, c("mahalanobis", "robust_mahalanobis", "scaled_euclidean", "euclidean"))

    distance.mat <- {
      if (distance.mat == "mahalanobis") {
        MatchIt::mahalanobis_dist(formula, dataset)
      }
      else if (distance.mat == "scaled_euclidean") {
        MatchIt::scaled_euclidean_dist(formula, dataset)
      }
      else if (distance.mat == "euclidean") {
        MatchIt::euclidean_dist(formula, dataset)
      }
      else if (distance.mat == "robust_mahalanobis") {
        MatchIt::robust_mahalanobis_dist(formula, dataset)
      }
    }
  }
  else {
    if (verbose) cat("Checking distance matrix...\n")

    if (inherits(distance.mat, "dist")) distance.mat <- as.matrix(distance.mat)

    okay <- is.matrix(distance.mat) && !any(distance.mat < 0)
    if (okay) {
      if (all(dim(distance.mat) == length(treat))) {
        if (!isSymmetric(unname(distance.mat))) okay <- FALSE
        else distance.mat <- distance.mat[treat == 1, treat == 0]
      }
      else if (nrow(distance.mat) != sum(treat == 1) || ncol(distance.mat) != sum(treat == 0)) {
        okay <- FALSE
      }
    }
    if (!okay) {
      customStop("'distance.mat' must be one of \"mahalanobis\", \"scaled\", or \"euclidean\" or a square, symmetric, N x N or N1 x N0 distance matrix when metric = \"dist\".",
                 "makeFrontier()")
    }
  }

  # if (verbose && metric != "custom") cat("Computing distance matrix...\n")
  #
  # if (metric == "mahal") {
  #   covs.mat <- get.covs.matrix(formula, dataset)
  #   distance.mat <- calculateMdist(covs.mat, treat)
  # }
  # else if (metric == "euclid") {
  #   covs.mat <- get.covs.matrix(formula, dataset)
  #   distance.mat <- calculateEdist(covs.mat, treat)
  # }
  # else if (metric == "custom") {
  #   if (is.null(distance.mat)) {
  #     customStop("'distance.mat' must be specified when metric = \"custom\".", "makeFrontier()")
  #   }
  #   if (!is.matrix(distance.mat) || !is.numeric(distance.mat) ||
  #       nrow(distance.mat) != sum(treat == 1) || ncol(distance.mat) != sum(treat == 0) ||
  #       anyNA(distance.mat)) {
  #     customStop("'distance.mat' must be an N1 x N0 numeric matrix when metric = \"custom\".", "makeFrontier()")
  #   }
  # }

  if (verbose) cat("Calculating frontier...\n")

  if (QOI == "FSATE") {
    frontier <- distToFrontierFSATE(distance.mat, treat, verbose)
  }
  else if (QOI == "FSATT") {
    frontier <- distToFrontierFSATT(distance.mat, treat, verbose)
  }

  if (verbose) cat("Done!\n")

  return(frontier)
}

# .DistFrontier <- function(treatment, dataset, formula, metric, QOI, distance.mat = NULL, verbose){
#
#   treat <- dataset[[treatment]]
#
#   if (verbose && metric != "custom") cat("Computing distance matrix...\n")
#
#   if (metric == "mahal") {
#     covs.mat <- get.covs.matrix(formula, dataset)
#     distance.mat <- calculateMdist(covs.mat, treat)
#   }
#   else if (metric == "euclid") {
#     covs.mat <- get.covs.matrix(formula, dataset)
#     distance.mat <- calculateEdist(covs.mat, treat)
#   }
#   else if (metric == "custom") {
#     if (is.null(distance.mat)) {
#       customStop("'distance.mat' must be specified when metric = \"custom\".", "makeFrontier()")
#     }
#     if (!is.matrix(distance.mat) || !is.numeric(distance.mat) ||
#         nrow(distance.mat) != sum(treat == 1) || ncol(distance.mat) != sum(treat == 0) ||
#         anyNA(distance.mat)) {
#       customStop("'distance.mat' must be an N1 x N0 numeric matrix when metric = \"custom\".", "makeFrontier()")
#     }
#   }
#
#   if (verbose) cat("Calculating frontier...\n")
#
#   if (QOI == "FSATE") {
#     frontier <- distToFrontierFSATE(distance.mat, treat, verbose)
#   }
#   else if (QOI == "FSATT") {
#     frontier <- distToFrontierFSATT(distance.mat, treat, verbose)
#   }
#
#   if (verbose) cat("Done!\n")
#
#   return(frontier)
# }
