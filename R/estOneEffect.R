estOneEffect <- function(formula, data, treatment, weights = NULL, QOI, subclass = NULL, id = NULL, alpha = .05) {

  data <- process_safe_for_model(data, formula)

  #For now, ignoring clustering by subclass and ID. Seems vcovHC is better for matching
  #w/ replacement.
  id <- NULL

  outcome <- model.response(model.frame(update(formula, . ~ 1), data = data))

  outcome.type <- if (length(unique(outcome)) == 2L) "binary" else "continuous"

  fit <- switch(outcome.type,
                "binary" = do.call("glm", list(formula, data = data,
                                               weights = weights,
                                               family = quasibinomial()),
                                   quote = TRUE),
                "continuous" = do.call("lm", list(formula, data = data,
                                                  weights = weights),
                                       quote = TRUE)
  )

  beta0 <- coef(fit)
  w <- weights(fit)

  if (QOI %in% c("SATT", "FSATT")) {
    if (!is.null(w)) w <- w[data[[treatment]] == 1]
    data <- data[data[[treatment]] == 1, , drop = FALSE]
  }

  estfun <- function(fit, beta = NULL) {
    if (!is.null(beta)) fit$coefficients <- beta

    data[[treatment]] <- 1

    p1 <- predict(fit, newdata = data, type = "response",
                  rankdeficient = "non-estim")

    Ep1 <- {
      if (is.null(w)) mean(p1)
      else weighted.mean(p1, w)
    }

    data[[treatment]] <- 0

    p0 <- predict(fit, newdata = data, type = "response",
                  rankdeficient = "non-estim")

    Ep0 <- {
      if (is.null(w)) mean(p0)
      else weighted.mean(p0, w)
    }

    Ep1 - Ep0
  }

  # Get estimate (difference in means)
  est <- {
    if (outcome.type == "continuous" &&
        length(beta0) == 2L &&
        all(names(beta0) %in% c("(Intercept)", treatment))) {
      beta0[treatment]
    }
    else estfun(fit)
  }

  if (is.null(alpha)) {
    return(setNames(est, "Estimate"))
  }

  # Get coefficient VCOV
  v <- {
    if (is.null(subclass) && is.null(id))
      safe_vcov(fit, vcov. = sandwich::vcovHC, type = "HC3")
    else if (!is.null(subclass))
      safe_vcov(fit, vcov. = sandwich::vcovCL, type = "HC1",
                cluster = list(subclass))
    else
      safe_vcov(fit, vcov. = sandwich::vcovCL, type = "HC1",
                cluster = list(subclass, id))
  }

  # Get estimate SE
  se <- NULL
  if (outcome.type == "continuous" &&
      abs(est - beta0[treatment]) < sqrt(.Machine$double.eps)) {
    se <- sqrt(v[treatment, treatment])
    if (!is.finite(se)) se <- NULL
  }

  if (is.null(se)) {
    #Get jacobian
    eps <- 1e-8
    J <- vapply(seq_along(beta0), function(i) {
      beta0[i] <- beta0[i] + eps
      (estfun(fit, beta0) - est) / eps
    }, numeric(1L))

    #Delta method SE
    se <- sqrt(rowSums(tcrossprod(J, v) * J))
  }

  #Get CI
  a <- c(alpha, 1 - alpha)
  fac <- qnorm(a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
  ci <- est + se * fac

  setNames(c(est, ci), c("Estimate", pct))
}

safe_vcov <- function(fit, vcov. = sandwich::vcovHC, type = "HC3", ...) {
  vcov.(fit, type = type, ...)
}