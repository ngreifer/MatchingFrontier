summary.matchFrontier <- function(object, N, Ndrop, ...) {

  if (missing(N) && missing(Ndrop)) {
    N <- NULL
  }
  else if (!missing(N) && !missing(Ndrop)) {
    customStop("only one of 'N' or 'Ndrop' may be specified.", "summary()")
  }
  else if (!missing(Ndrop)) {
    if (!is.numeric(Ndrop) ||
        any(Ndrop < 0) || any(Ndrop > object$n - 1)) {
      customStop(paste0("'Ndrop' must be between 0 and ", object$n - 1, "."), "summary()")
    }
    N <- object$n - Ndrop
  }
  else if (!missing(N)) {
    if (!is.numeric(N) ||
        any(N < 1) || any(N > object$n)) {
      customStop(paste0("'N' must be between 1 and ", object$n, "."), "summary()")
    }
  }

  if (is.null(N)) {
    startStat <- object$frontier$Ys[1]
    endStat <- object$frontier$Ys[length(object$frontier$Ys)]

    bestind <- which.min(object$frontier$Ys)
    bestStat <- object$frontier$Ys[bestind]

    treated.ind <- which(object$dataset[[object$treatment]] == 1)
    control.ind <- which(object$dataset[[object$treatment]] == 0)

    if (object$QOI %in% c("SATE", "FSATE")) {
      startNtreated <- length(treated.ind)
      startNcontrol <- length(control.ind)
      startN <- startNtreated + startNcontrol

      endNtreated <- sum(!treated.ind %in% unlist(object$frontier$drop.order))
      endNcontrol <- sum(!control.ind %in% unlist(object$frontier$drop.order))
      endN <- startN - object$frontier$Xs[length(object$frontier$Xs)]
      stopifnot(endNtreated + endNcontrol == endN)

      bestNtreated <- sum(!treated.ind %in% unlist(object$frontier$drop.order[seq_len(bestind)]))
      bestNcontrol <- sum(!control.ind %in% unlist(object$frontier$drop.order[seq_len(bestind)]))
      bestN <- startN - object$frontier$Xs[bestind]
    }
    else if (object$QOI %in% c("SATT")) {
      startNtreated <- length(treated.ind)
      startNcontrol <- length(control.ind)
      startN <- startNtreated + startNcontrol

      endNtreated <- startNtreated
      endNcontrol <- sum(!control.ind %in% unlist(object$frontier$drop.order))
      endN <- startN - object$frontier$Xs[length(object$frontier$Xs)]
      stopifnot(endNtreated + endNcontrol == endN)

      bestNtreated <- startNtreated
      bestNcontrol <- sum(!control.ind %in% unlist(object$frontier$drop.order[seq_len(bestind)]))
      bestN <- startN - object$frontier$Xs[bestind]
    }
    else if (object$QOI %in% c("FSATT")) {
      startNtreated <- length(treated.ind)
      startNcontrol <- NA
      startN <- startNtreated

      endNtreated <- sum(!treated.ind %in% unlist(object$frontier$drop.order))
      endNcontrol <- NA
      endN <- endNtreated

      bestNtreated <- sum(!treated.ind %in% unlist(object$frontier$drop.order[seq_len(bestind)]))
      bestNcontrol <- NA
      bestN <- bestNtreated
    }

    out <- list(
      Ntreated = c(start = startNtreated, end = endNtreated, best = bestNtreated),
      Ncontrol = c(start = startNcontrol, end = endNcontrol, best = bestNcontrol),
      N = c(start = startN, end = endN, best = bestN),
      Stat = c(start = startStat, end = endStat, best = bestStat),
      bestind = bestind,
      QOI = object$QOI,
      metric = object$metric)
  }
  else {
    ord <- order(N, decreasing = TRUE)
    N <- N[ord]

    selectedNtreated <- numeric(length(N))
    selectedNcontrol <- numeric(length(N))
    selectedN <- numeric(length(N))
    selectedStat <- numeric(length(N))

    Ns <- object$n - object$frontier$Xs

    treated.ind <- which(object$dataset[[object$treatment]] == 1)
    control.ind <- which(object$dataset[[object$treatment]] == 0)

    if (object$QOI %in% c("SATE", "FSATE", "SATT")) {
      startNtreated <- length(treated.ind)
      startNcontrol <- length(control.ind)
      startN <- startNtreated + startNcontrol
    }

    for (i in seq_along(N)) {
      ind <- which.min(Ns[Ns >= N[i]])
      N[i] <- Ns[ind]

      selectedStat[i] <- object$frontier$Ys[ind]

      if (object$QOI %in% c("SATE", "FSATE")) {
        selectedNtreated[i] <- sum(!treated.ind %in% unlist(object$frontier$drop.order[seq_len(ind)]))
        selectedNcontrol[i] <- sum(!control.ind %in% unlist(object$frontier$drop.order[seq_len(ind)]))
        selectedN[i] <- startN - object$frontier$Xs[ind]
      }
      else if (object$QOI %in% c("SATT")) {
        selectedNtreated[i] <- startNtreated
        selectedNcontrol[i] <- sum(!control.ind %in% unlist(object$frontier$drop.order[seq_len(ind)]))
        selectedN[i] <- startN - object$frontier$Xs[ind]
      }
      else if (object$QOI %in% c("FSATT")) {
        selectedNtreated[i] <- sum(!treated.ind %in% unlist(object$frontier$drop.order[seq_len(ind)]))
        selectedNcontrol[i] <- NA
        selectedN[i] <- selectedNtreated[i]
      }
    }

    names(selectedNtreated) <- if (missing(Ndrop)) paste("at N =", N) else paste("at Ndrop =", object$n - N)
    names(selectedNcontrol) <- names(selectedN) <- names(selectedStat) <- names(selectedNtreated)

    out <- list(
      Ntreated = selectedNtreated,
      Ncontrol = selectedNcontrol,
      N = selectedN,
      Stat = selectedStat,
      bestind = NULL,
      QOI = object$QOI,
      metric = object$metric)
  }

  class(out) <- "summary.matchFrontier"
  return(out)
}

print.summary.matchFrontier <- function(x, digits = 3, ...) {
  cat("Summary of matchFrontier object:\n\n")

  s <- matrix("", nrow = length(x[["N"]]), ncol = 4)


  s[!is.na(x[["Ntreated"]]), 1] <- as.character(as.integer(x[["Ntreated"]][!is.na(x[["Ntreated"]])]))
  s[!is.na(x[["Ncontrol"]]), 2] <- as.character(as.integer(x[["Ncontrol"]][!is.na(x[["Ncontrol"]])]))
  s[!is.na(x[["N"]]), 3] <- as.character(as.integer(x[["N"]][!is.na(x[["N"]])]))
  s[,4] <- format(x[["Stat"]], digits = digits, scientific = FALSE)

  s <- rbind(s, "")

  s[nrow(s), switch(x$QOI, "SATE" = 3, "FSATE" = 3, "SATT" = 2, "FSATT" = 1)] <- "^"

  colnames(s) <- c("N treated", "N control", "N total", "Statistic")
  rownames(s) <- c(firstup(names(x[["N"]])), "")

  print.data.frame(as.data.frame(s))

  invisible(x)
}
