makeMatchedData <- function(dataset, matched.to = NULL, drop.inds = NULL, weights = ".weights",
                            dup = FALSE, subclass = ".subclass", id = ".id") {
  if (is.null(matched.to)) {
    if (length(drop.inds) > 0) dataset <- dataset[-drop.inds,,drop=FALSE]

    dataset[[weights]] <- 1

    attr(dataset, "weights") <- weights
  }
  else if (dup) {
    if (length(drop.inds) != 0) {
      is.na(matched.to[drop.inds]) <- TRUE
    }

    s <- factor(c(which(!is.na(matched.to)), which(!is.na(matched.to))))
    levels(s) <- seq_len(nlevels(s))
    i <- c(which(!is.na(matched.to)), matched.to[!is.na(matched.to)])

    rn <- rownames(dataset)
    dataset <- dataset[i,]

    dataset <- setNames(cbind(rn[i], s, 1, dataset),
                        c(id, subclass, weights, names(dataset)))

    attr(dataset, "weights") <- weights
    attr(dataset, "subclass") <- subclass
    attr(dataset, "id") <- id
  }
  else {
    if (length(drop.inds) > 0) {
      is.na(matched.to)[drop.inds] <- TRUE
    }

    w <- as.numeric(!is.na(matched.to)) + tabulate(matched.to, nbins = length(matched.to))

    dataset[[weights]] <- w
    dataset <- dataset[w > 0,]

    attr(dataset, "weights") <- weights
  }
  return(dataset)
}
