#Remove units to lower energy distance between treated and full sample,
#between control and full sample, and between treated and control (Huling
#and Mak's "improved" energy distance)
energyToFrontierSATE <- function(distance.mat, treat.vec, verbose, keep.n.equal = FALSE) {
  #Energy distance formula:
  #2/n1n0 sum(D[t,c]) - 1/n1n1 sum(D[t,t]) - 1/n0n0 sum(D[c,c])

  N <- length(treat.vec)

  ind <- seq_along(treat.vec)
  treated.ind <- ind[treat.vec == 1]
  control.ind <- ind[treat.vec == 0]

  #Compute each unit's contribution to energy dist
  #Find unit that contributes maximally to energy distance
  #Remove that unit

  drop.order <- vector("list", N)
  Ys <- numeric(N)
  Xs <- integer(N)

  N1 <- length(treated.ind)
  N0 <- length(control.ind)

  if (verbose) {
    pb <- txtProgressBar(min = 0, max = N, style = 3)
  }

  Sdff <- sum(distance.mat)

  df1 <- distance.mat[,treated.ind, drop = FALSE]
  cSdf1 <- colSums(df1)
  Sdf1 <- sum(cSdf1)

  df0 <- distance.mat[,control.ind, drop = FALSE]
  cSdf0 <- colSums(df0)
  Sdf0 <- Sdff - Sdf1

  d10 <- df0[treated.ind,, drop = FALSE]
  rSd10 <- rowSums(d10)
  cSd10 <- colSums(d10)
  Sd10 <- sum(rSd10)

  d11 <- df1[treated.ind,,drop = FALSE]
  cSd11 <- colSums(d11)
  Sd11 <- sum(cSd11)

  d00 <- df0[control.ind,, drop = FALSE]
  cSd00 <- colSums(d00)
  Sd00 <- sum(cSd00)

  Ys[1] <- 2*Sd10/(N1*N0) - Sd11/(N1*N1) - Sd00/(N0*N0) +
    2*Sdf1/(N1*N) - Sd11/(N1*N1) - Sdff/(N*N) +
    2*Sdf0/(N0*N) - Sd00/(N0*N0) - Sdff/(N*N)

  min.edist <- Ys[1]

  #Difference between edist and edist with each unit
  treated.contributions <- (2/(N1*N0) - 2/((N1-1)*N0))*Sd10 +
    (-1/(N1^2) + 1/((N1-1)^2))*Sd11 +
    (2/((N1-1)*N0))*rSd10 -
    (2/((N1-1)^2))*cSd11 +
    (2/(N1*N) - 2/((N1-1)*N))*Sdf1 +
    (-1/(N1^2) + 1/((N1-1)^2))*Sd11 +
    (2/((N1-1)*N))*cSdf1 -
    (2/((N1-1)^2))*cSd11
  control.contributions <- (2/(N0*N1) - 2/((N0-1)*N1))*Sd10 +
    (-1/(N0^2) + 1/((N0-1)^2))*Sd00 +
    (2/((N0-1)*N1))*cSd10 -
    (2/((N0-1)^2))*cSd00 +
    (2/(N0*N) - 2/((N0-1)*N))*Sdf0 +
    (-1/(N0^2) + 1/((N0-1)^2))*Sd00 +
    (2/((N0-1)*N))*cSdf0 -
    (2/((N0-1)^2))*cSd00

  for (k in seq_len(N)[-1]) {

    #Find which units have the largest contribution to the energy distance
    if (keep.n.equal) {
      if (N0 > N1) {
        largest.contribution <- max(control.contributions)
        drop.treated <- integer(0L)
        drop.control <- which.min(abs(control.contributions - largest.contribution))
      }
      else {
        largest.contribution <- max(treated.contributions)
        drop.treated <- which.min(abs(treated.contributions - largest.contribution))
        drop.control <- integer(0L)
      }
    }
    else {
      largest.contribution <- max(max(treated.contributions), max(control.contributions))

      drop.treated <- which(abs(treated.contributions - largest.contribution) < 1e-9)
      drop.control <- which(abs(control.contributions - largest.contribution) < 1e-9)
    }

    treated.ind.to.drop <- treated.ind[drop.treated]
    control.ind.to.drop <- control.ind[drop.control]

    #Remove those units and decrease
    if (length(drop.treated) > 0) {
      treated.ind <- treated.ind[-drop.treated]
      N1 <- length(treated.ind)
    }
    if (length(drop.control) > 0) {
      control.ind <- control.ind[-drop.control]
      N0 <- length(control.ind)
    }

    #If removed units are last units, don't remove and stop
    if (N1-1 <= 0 || N0-1 <= 0) break

    if (verbose) {
      setTxtProgressBar(pb, N - (N0 + N1))
    }

    #Compute new edist with dropped units removed by subtracting contributions
    #of discarded units from remaining sum
    Sd10 <- Sd10 - sum(d10[drop.treated,]) -
      sum(d10[,drop.control]) -
      sum(d10[drop.treated,drop.control])

    if (length(drop.treated) > 0) {
      rSd10 <- rSd10[-drop.treated]
      cSd10 <- cSd10 - colSums(d10[drop.treated,, drop = FALSE])
      d10 <- d10[-drop.treated,, drop = FALSE]

      Sdf1 <- Sdf1 - sum(df1[,drop.treated])
      cSdf1 <- cSdf1[-drop.treated]
      df1 <- df1[, -drop.treated, drop = FALSE]

      Sd11 <- Sd11 - 2*sum(d11[-drop.treated, drop.treated]) -
        sum(d11[drop.treated, drop.treated])
      cSd11 <- cSd11[-drop.treated] - colSums(d11[drop.treated, -drop.treated, drop = FALSE])
      d11 <- d11[-drop.treated, -drop.treated, drop = FALSE]
    }
    if (length(drop.control) > 0) {
      rSd10 <- rSd10 - rowSums(d10[,drop.control, drop = FALSE])
      cSd10 <- cSd10[-drop.control]
      d10 <- d10[,-drop.control, drop = FALSE]

      Sdf0 <- Sdf0 - sum(df0[,drop.control])
      cSdf0 <- cSdf0[-drop.control]
      df0 <- df0[, -drop.control, drop = FALSE]

      Sd00 <- Sd00 - 2*sum(d00[-drop.control, drop.control]) -
        sum(d00[drop.control, drop.control])
      cSd00 <- cSd00[-drop.control] - colSums(d00[drop.control, -drop.control, drop = FALSE])
      d00 <- d00[-drop.control, -drop.control, drop = FALSE]
    }

    edist <- 2*Sd10/(N1*N0) - Sd11/(N1*N1) - Sd00/(N0*N0) +
      2*Sdf1/(N1*N) - Sd11/(N1*N1) - Sdff/(N*N) +
      2*Sdf0/(N0*N) - Sd00/(N0*N0) - Sdff/(N*N)

    #After passing sample size threshold of 90% of original N, stop if
    #new edist is larger than smallest edist
    if (edist < min.edist) min.edist <- edist
    else if ((N1+N0)/N < .9 && edist-min.edist > .2*(Ys[1]-min.edist)) break

    #Record new edist and units dropped
    Ys[k] <- edist
    drop.order[[k]] <- sort(c(treated.ind.to.drop, control.ind.to.drop))

    treated.contributions <- (2/(N1*N0) - 2/((N1-1)*N0))*Sd10 +
      (-1/(N1^2) + 1/((N1-1)^2))*Sd11 +
      (2/((N1-1)*N0))*rSd10 -
      (2/((N1-1)^2))*cSd11 +
      (2/(N1*N) - 2/((N1-1)*N))*Sdf1 +
      (-1/(N1^2) + 1/((N1-1)^2))*Sd11 +
      (2/((N1-1)*N))*cSdf1 -
      (2/((N1-1)^2))*cSd11
    control.contributions <- (2/(N0*N1) - 2/((N0-1)*N1))*Sd10 +
      (-1/(N0^2) + 1/((N0-1)^2))*Sd00 +
      (2/((N0-1)*N1))*cSd10 -
      (2/((N0-1)^2))*cSd00 +
      (2/(N0*N) - 2/((N0-1)*N))*Sdf0 +
      (-1/(N0^2) + 1/((N0-1)^2))*Sd00 +
      (2/((N0-1)*N))*cSdf0 -
      (2/((N0-1)^2))*cSd00
  }

  Ys <- Ys[seq_len(k-1)]
  drop.order <- drop.order[seq_len(k-1)]
  Xs <- cumsum(lengths(drop.order))

  if (verbose) {
    setTxtProgressBar(pb, N)
    close(pb)
  }

  return(list(drop.order = drop.order, Xs = Xs, Ys = Ys, Y.origin = Ys[1]))
}
