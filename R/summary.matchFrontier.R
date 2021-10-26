summary.matchFrontier <- function(object, ...) {

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

  class(out) <- "summary.matchFrontier"
  return(out)
}

print.summary.matchFrontier <- function(x, ...) {
  cat("Summary of matchFrontier object:\n\n")

  s <- matrix("", nrow = 3, ncol = 4)


  s[!is.na(x[["Ntreated"]]), 1] <- as.character(as.integer(x[["Ntreated"]][!is.na(x[["Ntreated"]])]))
  s[!is.na(x[["Ncontrol"]]), 2] <- as.character(as.integer(x[["Ncontrol"]][!is.na(x[["Ncontrol"]])]))
  s[!is.na(x[["N"]]), 3] <- as.character(as.integer(x[["N"]][!is.na(x[["N"]])]))
  s[,4] <- as.character(round(x[["Stat"]], 3))

  s <- rbind(s, "")

  s[4, switch(x$QOI, "SATE" = 3, "FSATE" = 3, "SATT" = 2, "FSATT" = 1)] <- "^"

  colnames(s) <- c("N treated", "N control", "N total", "Statistic")
  rownames(s) <- c("Start", "End", "Best", "")

  print.data.frame(as.data.frame(s))

  invisible(x)
}
