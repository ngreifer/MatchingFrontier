binsToFrontierFSATE <- function(strata, treat.vec, metric = "l1", verbose, keep.n.equal = FALSE) {

  N1 <- sum(treat.vec == 1)
  N0 <- sum(treat.vec == 0)
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
    pb <- txtProgressBar(min = 0, max = N, style = 3)
  }

  strataholder <- lapply(unique(strata), function(s) which(strata == s))

  treated.counts <- vapply(strataholder, function(x) {
    sum(treat.vec[x] == 1)
  }, numeric(1L))

  control.counts <- vapply(strataholder, function(x) {
    sum(treat.vec[x] == 0)
  }, numeric(1L))

  treated.props <- treated.counts/N1
  control.props <- control.counts/N0

  diffs <- treated.props - control.props
  # diffs <- get.diffs(strataholder, treat.vec, N1, N0)
  Ys[1] <- Lstat(diffs)

  min.Lstat <- Ys[1]

  k <- 1
  repeat {
    k <- k + 1

    if (!keep.n.equal || N1 > N0) {
      altered.treat.props <- treated.counts/(N1 - 1)
      L.drop.treated <- vapply(seq_along(strataholder), function(i) {
        atp <- altered.treat.props
        atp[i] <- (treated.counts[i] - 1)/(N1-1)
        diffs <- atp - control.props
        Lstat(diffs)
      }, numeric(1L))
      best.L.treated <- min(L.drop.treated)
      if (keep.n.equal) best.L.control <- Inf
    }

    if (!keep.n.equal || N1 <= N0) {
      altered.control.props <- control.counts/(N0 - 1)
      L.drop.control <- vapply(seq_along(strataholder), function(i) {
        acp <- altered.control.props
        acp[i] <- (control.counts[i] - 1)/(N0-1)
        diffs <- acp - treated.props
        Lstat(diffs)
      }, numeric(1L))
      if (keep.n.equal) best.L.treated <- Inf
      best.L.control <- min(L.drop.control)
    }

    if (best.L.treated < best.L.control) {
      drop.group <- 1
      drop.from <- which(L.drop.treated == best.L.treated)[1]
    }
    else {
      drop.group <- 0
      drop.from <- which(L.drop.control == best.L.control)[1]
    }

    dropped.element.ind <- which(treat.vec[strataholder[[drop.from]]] == drop.group)[1]

    drop <- strataholder[[drop.from]][dropped.element.ind]
    strataholder[[drop.from]] <- strataholder[[drop.from]][-dropped.element.ind]

    if (drop.group == 1) {
      N1 <- N1 - 1
      treated.counts[drop.from] <- treated.counts[drop.from] - 1
      treated.props <- treated.counts/N1
      new.Lstat <- best.L.treated
    }
    else {
      N0 <- N0 - 1
      control.counts[drop.from] <- control.counts[drop.from] - 1
      control.props <- control.counts/N0
      new.Lstat <- best.L.control
    }

    if (N1 == 0 || N0 == 0) break

    if (verbose) {
      setTxtProgressBar(pb, N - (N0 + N1))
    }

    # if (new.Lstat > Ys[k - 1]) break

    if (new.Lstat < min.Lstat) min.Lstat <- new.Lstat
    else if ((N1+N0) < .9*N && new.Lstat - min.Lstat > .2*(Ys[1] - min.Lstat)) break

    Ys[k] <- new.Lstat
    drop.order[[k]] <- drop

    if (new.Lstat == 0) break
  }

  keep <- c(1L, which(lengths(drop.order) > 0))
  drop.order[-keep] <- NULL
  Ys <- Ys[keep]

  Xs <- cumsum(lengths(drop.order))

  if (verbose) {
    setTxtProgressBar(pb, N)
    close(pb)
  }

  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, Y.origin = Ys[1]))
}
