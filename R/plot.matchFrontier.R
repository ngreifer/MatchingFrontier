plot.matchFrontier <- function(x, covs = NULL, stat = "std-diff", ...) {

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
    stat <- match_arg(stat, c("std-diff", "diff", "ks", "std-mean", "mean", "ks-target"))

    dataset <- as.data.frame(get.covs.matrix(data = x$dataset[covs]))
    cov.names <- names(dataset)

    if (stat == "std-diff") {
      sf <- switch(x$QOI,
                   "SATE" =, "FSATE" = sqrt(.5*(vapply(dataset[x$dataset[[x$treatment]] == 1,, drop = FALSE], var, numeric(1L)) +
                                                  vapply(dataset[x$dataset[[x$treatment]] == 0,, drop = FALSE], var, numeric(1L)))),
                   "SATT" =, "FSATT" = vapply(dataset[x$dataset[[x$treatment]] == 1,,drop = FALSE], sd, numeric(1L)))
      stat.fun <- function(v, t, d, w = NULL) {
        (w_m(d[[v]][d[[t]] == 1], w[d[[t]]==1]) -
           w_m(d[[v]][d[[t]] == 0], w[d[[t]]==0])) / sf[v]
      }

    }
    else if (stat == "diff") {
      stat.fun <- function(v, t, d, w = NULL) {
        (w_m(d[[v]][d[[t]] == 1], w[d[[t]]==1]) -
           w_m(d[[v]][d[[t]] == 0], w[d[[t]]==0]))
      }
    }
    else if (stat == "ks") {
      bin.vars <- vapply(dataset, function(V) length(unique(V)) <= 2, logical(1L))

      stat.fun <- function(v, t, d, w = NULL) {
        if (bin.vars[v]) {
          abs(w_m(d[[v]][d[[t]] == 1], w[d[[t]] == 1]) -
                w_m(d[[v]][d[[t]] == 0], w[d[[t]] == 0]))
        }
        else {
          w_ks(d[[v]], d[[t]], w)
        }
      }
    }
    else if (stat == "std-mean") {
      center <- switch(x$QOI,
                       "SATE" =, "FSATE" = colMeans(dataset),
                       "SATT" =, "FSATT" = colMeans(dataset[x$dataset[[x$treatment]] == 1,,drop = FALSE]))
      sf <- switch(x$QOI,
                   "SATE" =, "FSATE" = sqrt(.5*(vapply(dataset[x$dataset[[x$treatment]] == 1,,drop = FALSE], var, numeric(1L)) +
                                                  vapply(dataset[x$dataset[[x$treatment]] == 0,,drop = FALSE], var, numeric(1L)))),
                   "SATT" =, "FSATT" = vapply(dataset[x$dataset[[x$treatment]] == 1,,drop = FALSE], sd, numeric(1L)))

      stat.fun <- function(v, t, d, w = NULL) {
        (w_m(d[[v]], w) - center[v]) / sf[v]
      }
    }
    else if (stat == "mean") {
      stat.fun <- function(v, t, d, w = NULL) {
        w_m(d[[v]], w)
      }
    }
    else if (stat == "ks-target") {

      bin.vars <- vapply(dataset, function(V) length(unique(V)) <= 2, logical(1L))
      if (any(bin.vars)) {
          center <- switch(x$QOI,
                           "SATE" =, "FSATE" = colMeans(dataset[bin.vars]),
                           "SATT" =, "FSATT" = colMeans(dataset[x$dataset[[x$treatment]] == 1, bin.vars, drop = FALSE]))
      }

      stat.fun <- function(v, t, d, w = NULL) {
        if (bin.vars[v]) {
          abs(w_m(d[[v]], w) - center[v])
        }
        else {
          if (is.null(w)) w <- rep(1, nrow(d))
          switch(x$QOI,
                 "SATE" =, "FSATE" = w_ks(c(dataset[[v]], d[[v]]),
                                          treat = c(rep(1, nrow(dataset)), rep(0, nrow(d))),
                                          w = c(rep(1, nrow(dataset)), w)),
                 "SATT" =, "FSATT" = w_ks(c(dataset[[v]][dataset[[t]] == 1], d[[v]]),
                                          treat = c(rep(1, sum(dataset[[t]] == 1)), rep(0, nrow(d))),
                                          w = c(rep(1, sum(dataset[[t]] == 1)), w)))
        }
      }
    }

    Xs.long <- rep(x$frontier$Xs, ncol(dataset))
    dataset <- setNames(cbind(dataset, x$dataset[[x$treatment]]),
                        c(cov.names, x$treatment))

    cov.stat <- do.call("rbind", lapply(seq_along(x$frontier$Xs), function(i) {
      this.dat.inds.to.drop <- unlist(x$frontier$drop.order[seq_len(i)])

      m.data <- makeMatchedData(dataset,
                                matched.to = x$matched.to,
                                drop.inds = this.dat.inds.to.drop)

      vapply(cov.names, stat.fun, numeric(1L),
             t = x$treatment,
             d = m.data,
             w = m.data[[attr(m.data, "weights")]])
    }))

    cov.stat <- as.data.frame(cov.stat)

    cov.stat.long <- reshape(cov.stat, direction = "long",
                             varying = cov.names,
                             v.names = "Val",
                             timevar = "Covariate",
                             times = cov.names)

    #Unadjusted statistics
    cov.point.stat <- data.frame(Val = vapply(cov.names, stat.fun, numeric(1L),
                                              t = x$treatment,
                                              d = dataset),
                                 Covariate = cov.names
    )

    if (stat == "mean") {
      center <- switch(x$QOI,
                       "SATE" =, "FSATE" = colMeans(dataset[cov.names]),
                       "SATT" =, "FSATT" = colMeans(dataset[x$dataset[[x$treatment]] == 1, cov.names, drop = FALSE]))
      original.means <- data.frame(Covariate = cov.names, Val = center)

      p <- p + geom_hline(data = original.means,
                          aes(yintercept = Val, color = Covariate),
                          linetype = 2)
    }
    else {
      p <- p + geom_hline(yintercept = 0)
    }

    p <- p + plot_geom(data = cov.stat.long, mapping = aes(x = Xs.long, y = Val, color = Covariate)) +
      geom_point(data = cov.point.stat, mapping = aes(x = min(x$frontier$Xs), y = Val, color = Covariate))

    ylab <- switch(stat,
                   "std-diff" = "Standardized mean difference",
                   "diff" = "Mean difference",
                   "ks" = "Kolmogorov-Smirnov statistic",
                   "std-mean" = "Standardized mean",
                   "mean" = "Mean",
                   "ks-target" = "Sample-target Kolmogorov-\nSmirnov statistic")
  }

  p + scale_x_continuous(sec.axis = dup_axis(trans = ~ x$n - ., name = sub("dropped", "remaining", xlab),
                                             breaks = upper.breaks),
                         breaks = breaks) +
    labs(title = "Frontier Plot", x = xlab, y = ylab) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5))
}
