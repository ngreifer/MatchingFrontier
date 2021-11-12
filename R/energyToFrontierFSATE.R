#Remove units to lower between-group energy distance
energyToFrontierFSATE <- function(distance.mat, treat.vec, verbose) {
  #Energy distance formula:
  #2/n1n0 sum(D[t,c]) - 1/n1n1 sum(D[t,t]) - 1/n0n0 sum(D[c,c])

  N <- length(treat.vec)

  ind <- seq_along(treat.vec)
  treated.ind <- ind[treat.vec == 1]
  control.ind <- ind[treat.vec == 0]

  ind.in.vec <- integer(N)
  ind.in.vec[treated.ind] <- seq_along(treated.ind)
  ind.in.vec[control.ind] <- seq_along(control.ind)

  #Compute each unit's contribution to energy dist
  #Find unit that contributes maximally to energy distance
  #Remove that unit

  drop.order <- vector("list", N)
  Ys <- numeric(N)
  Xs <- integer(N)

  N1 <- length(treated.ind)
  N0 <- length(control.ind)

  if (verbose) {
    pb <- txtProgressBar(min = 1, max = N, style = 3)
  }

  d10 <- distance.mat[treated.ind, control.ind, drop = FALSE]
  Sd10 <- sum(d10)

  d11 <- distance.mat[treated.ind, treated.ind, drop = FALSE]
  Sd11 <- sum(d11)

  d00 <- distance.mat[control.ind, control.ind, drop = FALSE]
  Sd00 <- sum(d00)


  Ys[1] <- 2*Sd10/(N1*N0) - Sd11/(N1*N1) - Sd00/(N0*N0)

  min.edist <- Ys[1]

  treated.contributions <- (2/(N1*N0) - 2/((N1-1)*N0))*Sd10 +
    (-1/(N1^2) + 1/((N1-1)^2))*Sd11 +
    (2/((N1-1)*N0))*rowSums(d10) -
    (2/((N1-1)^2))*colSums(d11)
  control.contributions <- (2/(N0*N1) - 2/((N0-1)*N1))*Sd10 +
    (-1/(N0^2) + 1/((N0-1)^2))*Sd00 +
    (2/((N0-1)*N1))*colSums(d10) -
    (2/((N0-1)^2))*colSums(d00)


  if (verbose) setTxtProgressBar(pb, 1)

  for (k in seq_len(N)[-1]) {

    #Find which units have the largest contribution to the energy distance
    largest.contribution <- max(max(treated.contributions), max(control.contributions))

    drop.treated <- which(abs(treated.contributions - largest.contribution) < 1e-9)
    drop.control <- which(abs(control.contributions - largest.contribution) < 1e-9)

    treated.ind.to.drop <- treated.ind[drop.treated]
    control.ind.to.drop <- control.ind[drop.control]

    #Remove those units and decrease
    if (length(drop.treated) > 0) {
      treated.ind <- treated.ind[-drop.treated]
      N1 <- length(treated.ind)
    } else if (length(drop.control) > 0) {
      control.ind <- control.ind[-drop.control]
      N0 <- length(control.ind)
    }

    #If removed units are last units, don't remove and stop
    if (N1-1 <= 0 || N0-1 <= 0) break

    #Compute new edist with dropped units removed
    d10 <- distance.mat[treated.ind, control.ind, drop = FALSE]
    Sd10 <- sum(d10)
    if (length(drop.treated) > 0) {
      d11 <- distance.mat[treated.ind, treated.ind, drop = FALSE]
      Sd11 <- Sd11 - 2*sum(distance.mat[treated.ind.to.drop,treated.ind]) -
        sum(distance.mat[treated.ind.to.drop,treated.ind.to.drop])
    }
    if (length(drop.control) > 0) {
      d00 <- distance.mat[control.ind, control.ind, drop = FALSE]
      Sd00 <- Sd00 - 2*sum(distance.mat[control.ind.to.drop,control.ind]) -
        sum(distance.mat[control.ind.to.drop,control.ind.to.drop])
    }

    edist <- 2*Sd10/(N1*N0) - Sd11/(N1*N1) - Sd00/(N0*N0)

    #After passing sample size threhsold of 90% of original N, stop if
    #new edist is larger than smallest edist
    if (edist < min.edist) min.edist <- edist
    else if ((N1+N0) < .9*N && edist-min.edist > .2*(Ys[1]-min.edist)) break

    #Record new edist and units dropped
    Ys[k] <- edist
    drop.order[[k]] <- sort(c(treated.ind.to.drop, control.ind.to.drop))

    #Compute unit contributions after having dropped units
    treated.contributions <- (2/(N1*N0) - 2/((N1-1)*N0))*Sd10 +
      (-1/(N1^2) + 1/((N1-1)^2))*Sd11 +
      (2/((N1-1)*N0))*rowSums(d10) -
      (2/((N1-1)^2))*colSums(d11)
    control.contributions <- (2/(N0*N1) - 2/((N0-1)*N1))*Sd10 +
      (-1/(N0^2) + 1/((N0-1)^2))*Sd00 +
      (2/((N0-1)*N1))*colSums(d10) -
      (2/((N0-1)^2))*colSums(d00)

    if (verbose) {
      setTxtProgressBar(pb, k)
    }
  }

  Ys <- Ys[seq_len(k-1)]
  drop.order <- drop.order[seq_len(k-1)]
  Xs <- cumsum(lengths(drop.order))

  if (verbose) {
    setTxtProgressBar(pb, N)
    close(pb)
  }

  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, distance.mat = distance.mat))
}
