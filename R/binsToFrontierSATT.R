binsToFrontierSATT <- function(strata, treat.vec, metric = "l1", verbose, keep.n.equal = FALSE){

  N1 <- sum(treat.vec == 1)
  N0 <- N0_ <- sum(treat.vec == 0)
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
    pb <- txtProgressBar(min = 0, max = N0, style = 3)
  }

  strataholder <- lapply(unique(strata), function(s) which(strata == s))
  diffs <- get.diffs(strataholder, treat.vec, N1, N0)
  Ys[1] <- Lstat(diffs)

  min.Lstat <- Ys[1]

  k <- 1
  repeat {
    k <- k + 1

    drop.from <- which.max(diffs)
    dropped.element.ind <- which(treat.vec[strataholder[[drop.from]]] == 0)[1]

    drop <- strataholder[[drop.from]][dropped.element.ind]
    strataholder[[drop.from]] <- strataholder[[drop.from]][-dropped.element.ind]

    N0 <- N0 - 1

    if (N0 == 0) break

    if (verbose) {
      setTxtProgressBar(pb, N0_ - N0)
    }

    diffs <- get.diffs(strataholder, treat.vec, N1, N0)
    new.Lstat <- Lstat(diffs)

    if (new.Lstat < min.Lstat) min.Lstat <- new.Lstat
    else if (N0 < .9*(N-N1) && new.Lstat-min.Lstat > .2*(Ys[1]-min.Lstat)) break

    Ys[k] <- new.Lstat
    drop.order[[k]] <- drop
  }

  keep <- c(1L, which(lengths(drop.order) > 0))
  drop.order[-keep] <- NULL
  Ys <- Ys[keep]

  Xs <- cumsum(lengths(drop.order))

  if (verbose) {
    setTxtProgressBar(pb, N0_)
    close(pb)
  }

  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, Y.origin = Ys[1]))
}
