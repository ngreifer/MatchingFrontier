#Find matches for treated units, then discard treated units with greatest distance
#matched.to is NA for control units
distToFrontierFSATT <- function(distance.mat, treat.vec, verbose) {

    treated.ind <- which(treat.vec == 1)
    control.ind <- which(treat.vec == 0)

    N <- length(treat.vec)
    N1 <- length(treated.ind)

    if (verbose) {
        pb <- txtProgressBar(min = 1, max = N1, style = 3)
    }

    #Find closest matches to treated units (rows)
    row.mins.inds <- apply(distance.mat, 1, function(x) which.min(x))

    #Assign index of who is matched to who
    matched.to <- control.ind[row.mins.inds]

    #Extract distances to closest matches
    min.distances <- distance.mat[cbind(seq_len(nrow(distance.mat)), row.mins.inds)]
    min.distances.full <- rep(NA, N)
    min.distances.full[treated.ind] <- min.distances

    inds <- seq_len(N1)

    #Order (rank) distances to closest matches
    ranks <- rank(min.distances, ties.method = "min", na.last = NA)

    #Create list of drop order, starting with dropping no one
    #Entries will be empty when there are ties (which are grouped together)
    drop.order <- c(list(integer(0)), lapply(rev(inds), function(i) treated.ind[ranks == i]))

    #Empty entries to remove (except no drop)
    empty <- c(which(lengths(drop.order) == 0)[-1], N1 + 1)

    #Sort distances to closest matches from smallest to largest
    sorted.min.distances <- c(0, rev(min.distances.full[unlist(drop.order)]))

    #Xs: how many have been removed at each step
    Xs <- c(0, inds)[-empty]

    #Ys: mean of remaining distances at each step
    Ys <- rev(cumsum(sorted.min.distances))[-empty] / (N1 - Xs)

    drop.order[empty] <- NULL

    matched.to.full <- rep(NA, N)
    matched.to.full[treated.ind] <- matched.to

    if (verbose) {
        setTxtProgressBar(pb, N1)
        close(pb)
    }

    # Checks to confirm monotonically decreasing. Since
    # that's theoretically impossible, if the condition is
    # met, there's a serious bug somewhere in the code.
    if (any(diff(Ys) > 0 )){
        stop('Something is very wrong. Email ngreifer@iq.harvard.edu.')
    }
    return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, matched.to = matched.to.full, distance.mat = distance.mat))
}

