frontier_to_matchit <- function(frontier.object, N, Ndrop) {

  if (!inherits(frontier.object, "matchFrontier")) {
    customStop("'frontier.object' must be a matchFrontier object, the output of a call to makeFrontier().")
  }
  if (missing(N) && missing(Ndrop)) {
    customStop("one of 'N' or 'Ndrop' must be specified.")
  }
  else if (!missing(N) && !missing(Ndrop)) {
    customStop("only one of 'N' or 'Ndrop' may be specified.")
  }
  else if (!missing(Ndrop)) {
    if (!is.numeric(Ndrop) || length(Ndrop) != 1 ||
        Ndrop < 0 || Ndrop > frontier.object$n - 1) {
      customStop(paste0("'Ndrop' must be a single number between 0 and ", frontier.object$n - 1, "."))
    }
    N <- frontier.object$n - Ndrop
  }
  else if (!missing(N)) {
    if (!is.numeric(N) || length(N) != 1 ||
        N < 1 || N > frontier.object$n) {
      customStop(paste0("'N' must be a single number between 1 and ", frontier.object$n, "."))
    }
  }

  Ns <- frontier.object$n - frontier.object$frontier$Xs
  ind <- which.min(Ns[Ns >= N])
  drop.inds <- unlist(frontier.object$frontier$drop.order[seq_len(ind)])

  estimand <- switch(frontier.object$QOI,
                     "SATE" = "ATE",
                     "FSATE" = "ATE",
                     "SATT" = "ATT",
                     "FSATT" = "ATT")

  lab <- rownames(frontier.object$dataset)
  treat <- setNames(frontier.object$dataset[[frontier.object$treatment]], lab)

  weights <- setNames(rep(0, nrow(frontier.object$dataset)), lab)

  if (inherits(frontier.object, "distFrontier")) {
    is.na(frontier.object$matched.to)[drop.inds] <- TRUE

    match.matrix <- matrix(lab[frontier.object$matched.to], ncol = 1,
                           dimnames = list(lab, NULL))

    weights[rownames(match.matrix)[!is.na(match.matrix[,1])]] <- 1
    for (i in unique(match.matrix[!is.na(match.matrix[,1]),1])) weights[i] <- weights[i] + sum(!is.na(match.matrix[,1]) & match.matrix[,1] == i)
  }
  else {
    match.matrix <- NULL
    if (frontier.object$QOI %in% c("SATE", "FSATE", "FSATT")) weights[-drop.inds] <- 1
    else if (frontier.object$QOI %in% c("SATT")) {
      weights[treat == 1 | (treat == 0 & !seq_along(treat) %in% drop.inds)] <- 1
    }
  }

  #Compute sample sizes
  ESS <- function(w) sum(w)^2/sum(w^2)
  nn <- matrix(c(sum(treat==0), sum(treat==1),
                sum(treat==0), sum(treat==1),
                ESS(weights[treat==0]), ESS(weights[treat==1]),
                sum(treat==0 & weights > 0), sum(treat==1 & weights > 0),
                sum(treat==0 & weights==0), sum(treat==1 & weights==0),
                0, 0), byrow = TRUE, ncol=2, nrow=6,
              dimnames = list(c("All (ESS)", "All", "Matched (ESS)","Matched", "Unmatched","Discarded"),
                              c("Control", "Treated")))

  method <- "matchingFrontier"
  attr(method, "method") <- paste0("matching frontier (minimizing the ", metric2info(frontier.object$metric), ")")
  info <- list(method = method,
               replace = if (inherits(frontier.object, "distFrontier")) TRUE else NULL,
               mahalanobis = FALSE,
               distance_is_matrix = frontier.object$metric == "custom")
  if (inherits(frontier.object, "distFrontier")) {
    info$distance <- "user"
    attr(info$distance, "custom") <- switch(frontier.object$metric,
                                            "mahal" = "Mahalanobis",
                                            "euclid" = "Euclidean",
                                            NULL)
  }

  out <- list(
    match.matrix = match.matrix,
    weights = weights,
    X = frontier.object$dataset[frontier.object$match.on],
    call = frontier.object$call,
    info = info,
    estimand = estimand,
    formula = reformulate(frontier.object$match.on, frontier.object$treatment),
    treat = treat,
    nn = nn)

  class(out) <- "matchit"
  return(out)

}
