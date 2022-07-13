makeFrontier <- function(x, ...) {
  UseMethod("makeFrontier")
}
makeFrontier.data.frame <- function(x, treatment, match.on, QOI = 'FSATT',
                                    metric = 'dist', breaks = NULL,
                                    distance.mat = NULL, keep.n.equal = FALSE,
                                    verbose = TRUE, ...){

  call <- match.call()
  call[[1]] <- quote(makeFrontier)

  # Check the frontier arguments
  processed_metric <- processMetric(metric, breaks, distance.mat)

  checkArgs(QOI, metric = processed_metric, ...)

  if (!missing(match.on) && inherits(match.on, "formula")) {
    match.on <- all.vars(delete.response(terms(match.on)))
  }

  # Check data
  checkDat(x, treatment, match.on)

  x[[treatment]] <- binarize(x[[treatment]])

  formula <- reformulate(match.on, treatment)

  makeFrontier_internal(x, treatment, match.on, formula, QOI,
                        processed_metric, call, verbose)
}

makeFrontier.formula <- function(formula, data, QOI = 'FSATT',
                                 metric = 'dist', breaks = NULL,
                                 distance.mat = NULL, keep.n.equal = FALSE,
                                 verbose = TRUE, ...){

  call <- match.call()
  call[[1]] <- quote(makeFrontier)

  # Check the frontier arguments
  processed_metric <- processMetric(metric, breaks, distance.mat)

  checkArgs(QOI, metric = processed_metric, ...)

  treatment <- all.vars(formula)[attr(terms(formula), "response")]

  if (missing(treatment) || length(treatment) == 0){
    customStop("the treatment must be supplied.", 'makeFrontier()')
  }

  match.on <- setdiff(all.vars(formula), treatment)

  checkDat(data, treatment, match.on)

  data[[treatment]] <- binarize(data[[treatment]])

  makeFrontier_internal(data, treatment, match.on, formula, QOI,
                        processed_metric, call, keep.n.equal, verbose)
}

makeFrontier_internal <- function(dataset, treatment, match.on, formula, QOI = 'FSATT',
                                  metric, call = NULL, keep.n.equal = FALSE, verbose = FALSE) {


  frontier <- switch(attr(metric, "type"),
                     "dist" = DistFrontier(treatment = treatment,
                                           dataset = dataset,
                                           formula = formula,
                                           metric = metric,
                                           QOI = QOI,
                                           verbose = verbose),
                     "bin" = BinFrontier(treatment = treatment,
                                         dataset = dataset,
                                         formula = formula,
                                         metric = metric,
                                         QOI = QOI,
                                         match.on = match.on,
                                         keep.n.equal = keep.n.equal,
                                         verbose = verbose),
                     "energy" = EnergyFrontier(treatment = treatment,
                                               dataset = dataset,
                                               formula = formula,
                                               metric = metric,
                                               QOI = QOI,
                                               keep.n.equal = keep.n.equal,
                                               verbose = verbose)
  )

  matched.to <- frontier$matched.to
  frontier$matched.to <- NULL

  out <- list(
    frontier = frontier,
    treatment = treatment,
    QOI = QOI,
    metric = as.character(metric),
    dataset = dataset,
    match.on = match.on,
    matched.to = matched.to,
    call = call,
    n = switch(QOI,
               "SATE" = nrow(dataset),
               "FSATE" = nrow(dataset),
               "SATT" = sum(dataset[[treatment]] == 0),
               "FSATT" = sum(dataset[[treatment]] == 1))
  )

  class(out) <- c(paste0(attr(metric, "type"), "Frontier"), "matchFrontier")

  return(out)
}

print.matchFrontier <- function(x, ...){
  cat("A matchFrontier object\n")
  cat(paste0("- quantity of interest: ", x$QOI, "\n"))
  cat(paste0("- imbalance metric: ", metric2info(x$metric), "\n"))
  cat(paste0("- treatment: ", x$treatment, "\n"))
  cat(paste0("- covariates: ", paste0(x$match.on, collapse = ", "), "\n"))
  cat(paste0("- number of points: ", length(x$frontier$Xs), "\n"))

  invisible(x)
}
