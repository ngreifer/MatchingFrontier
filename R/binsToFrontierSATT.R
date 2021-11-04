binsToFrontierSATT <- function(strataholder, treat.vec, metric = "L1"){

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

    diffs <- get.diffs(strataholder, treat.vec, N1, N0)
    new.Lstat <- Lstat(diffs)


    if (new.Lstat < min.Lstat) min.Lstat <- new.Lstat
    else if (N0 < .9*(N-N1) && new.Lstat-min.Lstat > .2*(Ys[1]-min.Lstat)) break

    Ys[k] <- new.Lstat
    drop.order[[k]] <- drop

  }

  empty <- which(lengths(drop.order) == 0)
  drop.order[empty] <- NULL
  Ys <- Ys[-empty]

  Xs <- cumsum(lengths(drop.order))
  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys))
}
