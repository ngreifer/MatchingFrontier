binsToFrontierFSATE <- function(strataholder, treat.vec, metric = "l1", verbose) {

  N1 <- sum(treat.vec[unlist(strataholder)] == 1)
  N0 <- sum(treat.vec[unlist(strataholder)] == 0)
  N <- length(treat.vec)

  drop.order <- vector("list", N0)
  Ys <- numeric(N0)

  if (startsWith(metric, "l1")) {
    Lstat <- function(diffs) .5*sum(abs(diffs))
  }
  else if (startsWith(metric, "l2")) {
    Lstat <- function(diffs) .5*sqrt(sum(diffs^2))
  }

  if (verbose) {
    pb <- txtProgressBar(min = 1, max = N, style = 3)
  }

  diffs <- get.diffs(strataholder, treat.vec, N1, N0)
  Ys[1] <- Lstat(diffs)

  min.Lstat <- Ys[1]

  if (verbose) setTxtProgressBar(pb, 1)

  k <- 1
  repeat {
    k <- k + 1

    #Drop unit from group with higher proportion in stratum with
    #largest |diff|
    drop.from <- which.max(abs(diffs))

    if (diffs[drop.from] < 0) drop.group <- 1
    else drop.group <- 0

    dropped.element.ind <- which(treat.vec[strataholder[[drop.from]]] == drop.group)[1]

    drop <- strataholder[[drop.from]][dropped.element.ind]
    strataholder[[drop.from]] <- strataholder[[drop.from]][-dropped.element.ind]

    if (drop.group == 1) N1 <- N1 - 1
    else N0 <- N0 - 1

    if (N1 == 0 || N0 == 0) break

    diffs <- get.diffs(strataholder, treat.vec, N1, N0)
    new.Lstat <- Lstat(diffs)

    # if (new.Lstat > Ys[k - 1]) break

    if (new.Lstat < min.Lstat) min.Lstat <- new.Lstat
    else if ((N1+N0) < .9*N && new.Lstat - min.Lstat > .2*(Ys[1] - min.Lstat)) break

    Ys[k] <- new.Lstat
    drop.order[[k]] <- drop

    if (verbose) {
      setTxtProgressBar(pb, k)
    }

  }

  empty <- which(lengths(drop.order) == 0)
  drop.order[empty] <- NULL
  Ys <- Ys[-empty]

  Xs <- cumsum(lengths(drop.order))

  if (verbose) {
    setTxtProgressBar(pb, N)
    close(pb)
  }

  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys))
}
