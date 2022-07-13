#Find matches for all units, then discard units with greatest distance
distToFrontierFSATE <- function(distance.mat, treat.vec, verbose){

  N <- length(treat.vec)

  treated.ind <- which(treat.vec == 1)
  control.ind <- which(treat.vec == 0)

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = N, style = 3)
  }

  #Find closest matches to treated units (rows) and control units (cols)
  row.mins.inds <- apply(distance.mat, 1, function(x) which.min(x))
  col.mins.inds <- apply(distance.mat, 2, function(x) which.min(x))

  #Assign index of who is matched to who
  matched.to <- integer(N)
  matched.to[treated.ind] <- control.ind[row.mins.inds]
  matched.to[control.ind] <- treated.ind[col.mins.inds]

  #Extract distances to closest matches
  min.distances <- integer(N)
  min.distances[treated.ind] <- distance.mat[cbind(seq_len(nrow(distance.mat)), row.mins.inds)]
  min.distances[control.ind] <- distance.mat[cbind(col.mins.inds, seq_len(ncol(distance.mat)))]

  inds <- seq_along(treat.vec)

  #Order (rank) distances to closest matches
  ranks <- rank(min.distances, ties.method = "min", na.last = NA)

  #Create list of drop order, starting with dropping no one
  #Entries will be empty when there are ties (which are grouped together)
  drop.order <- c(list(integer(0)), lapply(c(rev(inds)), function(i) which(ranks == i)))

  #Empty entries to remove (except no drop)
  empty <- c(which(lengths(drop.order) == 0)[-1], N + 1)

  #Sort distances to closest matches from smallest to largest
  sorted.min.distances <- c(0, rev(min.distances[unlist(drop.order)]))

  #Xs: how many have been removed at each step
  Xs <- c(0, inds)[-empty]

  #Ys: mean of remaining distances at each step
  Ys <- rev(cumsum(sorted.min.distances))[-empty] / (N - Xs)

  drop.order[empty] <- NULL

  if (verbose) {
    setTxtProgressBar(pb, N)
    close(pb)
  }

  # Checks to confirm monotonically decreasing. Since
  # that's theoretically impossible, if the condition is
  # met, there's a serious bug somewhere in the code.
  if(any(diff(Ys) > 0 )){
    stop('Something is very wrong. Email ngreifer@iq.harvard.edu.')
  }
  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, matched.to = matched.to))
}

