plot.matchFrontier <- function(x, covs = NULL, diffs = TRUE, std = TRUE, ...) {

  xlab <- switch(x$QOI,
                 "SATE" = "Number of units dropped",
                 "FSATE" = "Number of units dropped",
                 "SATT" = "Number of control units dropped",
                 "FSATT" = "Number of treated units dropped")

  ylab <- firstup(metric2info(x$metric))

  breaks <- scales::breaks_extended()(x$frontier$Xs)
  upper.breaks <- x$n - breaks

  plot_geom <- if (length(x$frontier$Xs) == 1) geom_point else geom_line

  p <- ggplot()

  if (is.null(covs)) {
    p <- p + plot_geom(data = NULL, mapping = aes(x = x$frontier$Xs, y = x$frontier$Ys), ...)
  }
  else {
    dataset <- get.covs.matrix(data = x$dataset[covs])
    cov.names <- colnames(dataset)
    if (x$QOI %in% c("SATE", "FSATE")) {
      center <- colMeans(dataset)
    }
    else if (x$QOI %in% c("SATT", "FSATT")) {
      center <- colMeans(dataset[x$dataset[[x$treatment]] == 1,,drop = FALSE])
    }

    if (diffs || std) {
      p <- p + geom_hline(yintercept = 0)
    }
    else {
      original.means <- data.frame(Covariate = cov.names, Val = center)
      p <- p + geom_hline(data = original.means, aes(yintercept = Val, color = Covariate),
                          linetype = 2)
    }

    if (std) {
      if (x$QOI %in% c("SATE", "FSATE")) {
        sds <- apply(dataset, 2, sd)
      }
      else if (x$QOI %in% c("SATT", "FSATT")) {
        sds <- apply(dataset[x$dataset[[x$treatment]] == 1,,drop = FALSE], 2, sd)
      }
      dataset <- scale(dataset, center = center, scale = sds)
    }
    Xs.long <- rep(x$frontier$Xs, ncol(dataset))
    dataset <- setNames(cbind(as.data.frame(dataset), x$dataset[[x$treatment]]),
                        c(cov.names, x$treatment))

    if (diffs) {
      cov.stat <- do.call("rbind", lapply(seq_along(x$frontier$Xs), function(i) {
        this.dat.inds.to.drop <- unlist(x$frontier$drop.order[seq_len(i)])

        matched.dataset <- makeMatchedData(dataset,
                                           matched.to = x$matched.to,
                                           drop.inds = this.dat.inds.to.drop)

        w <- matched.dataset[[attr(matched.dataset, "weights")]]
        treat <- matched.dataset[[x$treatment]]

        vapply(matched.dataset[treat == 1, cov.names, drop = FALSE], weighted.mean,
               numeric(1L), w = w[treat == 1]) -
          vapply(matched.dataset[treat == 0, cov.names, drop = FALSE], weighted.mean,
                 numeric(1L), w = w[treat == 0])
      }))

      cov.point.stat <- data.frame(
        Val = colMeans(dataset[dataset[[x$treatment]] == 1, cov.names]) -
        colMeans(dataset[dataset[[x$treatment]] == 0, cov.names]),
        Covariate = cov.names
      )
    }
    else {
      cov.stat <- do.call("rbind", lapply(seq_along(x$frontier$Xs), function(i) {
        this.dat.inds.to.drop <- unlist(x$frontier$drop.order[seq_len(i)])

        matched.dataset <- makeMatchedData(dataset,
                                           matched.to = x$matched.to,
                                           drop.inds = this.dat.inds.to.drop)

        w <- matched.dataset[[attr(matched.dataset, "weights")]]

        vapply(matched.dataset[,cov.names, drop = FALSE], weighted.mean,
               numeric(1L), w = w)
      }))

      cov.point.stat <- NULL
    }

    cov.stat <<- as.data.frame(cov.stat)
    cov.stat <- as.data.frame(cov.stat)

    cov.stat.long <- reshape(cov.stat, direction = "long",
                             varying = cov.names,
                             v.names = "Val",
                             timevar = "Covariate",
                             times = cov.names)

    p <- p + plot_geom(data = cov.stat.long, mapping = aes(x = Xs.long, y = Val, color = Covariate))
    if (!is.null(cov.point.stat)) {
      p <- p + geom_point(data = cov.point.stat, mapping = aes(x = min(x$frontier$Xs), y = Val, color = Covariate))
    }

    if (diffs && std) ylab <- "Standardized mean difference"
    else if (diffs && !std) ylab <- "Mean difference"
    else if (!diffs && std) ylab <- "Standardized mean"
    else ylab <- "Mean"
  }

  p + scale_x_continuous(sec.axis = dup_axis(trans = ~ x$n - ., name = sub("dropped", "remaining", xlab),
                                             breaks = upper.breaks),
                         breaks = breaks) +
    labs(title = "Frontier Plot", x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))
}
