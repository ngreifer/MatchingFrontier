getBins <- function(dataset, match.on, breaks = NULL){

  # Remove dropped covs
  dataset <- dataset[match.on]
  for (i in seq_len(ncol(dataset))) {
    if (is.character(dataset[[i]])) dataset[[i]] <- factor(dataset[[i]])
    if (is.factor(dataset[[i]])) dataset[[i]] <- as.numeric(dataset[[i]])
  }

  #Process breaks; adapted from MatchIt:::cem_matchit()
  if (is.null(breaks)) breaks <- "sturges"

  if (!is.list(breaks)) {
    breaks <- setNames(as.list(rep(breaks, length(match.on))), match.on)
  }

  if (is.null(names(breaks))) customStop("'breaks' must be a named list of binning values with an element for each numeric variable.", "makeFrontier()")
  bad.names <- setdiff(names(breaks), match.on)
  nb <- length(bad.names)
  if (nb > 0) {
    breaks[bad.names] <- NULL
  }

  bad.breaks <- setNames(rep(FALSE, length(breaks)), names(breaks))
  for (i in names(breaks)) {
    if (length(breaks[[i]]) == 0) {
      breaks[[i]] <- "sturges"
    }
    else if (length(breaks[[i]]) == 1) {
      if (is.character(breaks[[i]])) {

        bad.breaks[i] <- is.na(pmatch(breaks[[i]], c("sturges", "fd", "scott")))
      }
      else if (is.numeric(breaks[[i]])) {
        if      (!is.finite(breaks[[i]]) || breaks[[i]] < 0) bad.breaks[i] <- TRUE
        if      (breaks[[i]] %in% 0:1) match.on <- setdiff(match.on, i) #Will not be binned
      }
    }
    else {
      bad.breaks[i] <- TRUE
    }
  }
  if (any(bad.breaks)) {
    customStop("each entry in the list supplied to 'breaks' must be a string containing the name of an allowable binning method or a single number corresponding to the number of bins",
               "makeFrontier()")
  }

  #Create bins for numeric variables
  for (i in match.on) {
    if (!i %in% names(breaks)) {
      if (is.character(dataset[[i]]) || is.factor(dataset[[i]])) next
      else breaks[[i]] <- "sturges"
    }

    if (is.character(breaks[[i]])) {
      breaks[[i]] <- match_arg(tolower(breaks[[i]]), c("sturges", "fd", "scott"))
      breaks[[i]] <- switch(breaks[[i]],
                            sturges = nclass.Sturges(dataset[[i]]),
                            fd = nclass.FD(dataset[[i]]),
                            scott = nclass.scott(dataset[[i]]))
      #Breaks is now a single number
    }

    #breaks is number of bins
    bins <- seq(min(dataset[[i]]), max(dataset[[i]]), length = breaks[[i]] + 1)
    bins[c(1, length(bins))] <- c(-Inf, Inf)

    dataset[[i]] <- findInterval(dataset[[i]], bins)
  }

  # Calculate strata
  strata <- factor(do.call("paste", c(dataset, list(sep = "|"))))

  strataholder <- lapply(levels(strata), function(u) which(strata==u))

  return(strataholder)
}

get.diffs <- function(strataholder, treat.vec, num.treated, num.control){
  #Diffs positive when control proportion larger than treated proportion
  unlist(lapply(strataholder, function(x)
    sum(treat.vec[x] == 0) / num.control - sum(treat.vec[x] == 1) / num.treated))
}

getBinsAtMedian <- function(dataset, match.on, treat.vec, metric){
  if (startsWith(metric, "L1")) {
    Lstat <- function(diffs) .5*sum(abs(diffs))
  }
  else if (startsWith(metric, "L2")) {
    Lstat <- function(diffs) .5*sqrt(sum(diffs^2))
  }

  num.treated <- sum(treat.vec == 1)
  num.control <- sum(treat.vec == 0)

  n.coarsenings <- 1001

  dataset <- dataset[match.on]
  for (i in seq_len(ncol(dataset))) {
    if (is.character(dataset[[i]])) dataset[[i]] <- factor(dataset[[i]])
    dataset[[i]] <- as.numeric(dataset[[i]])
  }

  nunique.covs <- vapply(match.on, function(i) length(unique(dataset[[i]])), integer(1L))

  strataholders <- lapply(seq_len(n.coarsenings), function(i) {
    breaks <- setNames(lapply(nunique.covs, function(nu) {
      sample(seq(min(2, nu), min(12, nu)), 1)
    }), match.on)

    getBins(dataset, match.on, breaks)
  })

  Lstats <- vapply(strataholders, function(s) {
    Lstat(get.diffs(s, treat.vec, num.treated, num.control))
  }, numeric(1L))

  Lstat.med <- sort(Lstats, partial = ceiling(n.coarsenings/2))[ceiling(n.coarsenings/2)]

  return(strataholders[[which(Lstats == Lstat.med)[1]]])

}
