checkDat <- function(dataset, treatment, match.on){

    if (missing(treatment) || length(treatment) == 0 || !is.character(treatment)) {
        customStop("a treatment variable must be supplied.", 'makeFrontier()')
    }

    if (missing(match.on) || length(match.on) == 0 || !is.character(match.on)) {
        customStop("matching variables must be supplied.", 'makeFrontier()')
    }

    if (missing(dataset) || length(dataset) == 0) {
        customStop("a dataset must be supplied.", 'makeFrontier()')
    }

    # Make sure user isn't trying to match on the treatment or the outcome
    if (treatment %in% match.on) {
        customStop("the treatment must not be in the matching variables.", 'makeFrontier()')
    }

    if (!all(match.on %in% names(dataset))) {
        customStop("all matching variables must be present in 'dataset'.", 'makeFrontier()')
    }

    if (!treatment %in% names(dataset)) {
        customStop("the treatment variable must be present in 'dataset'.", 'makeFrontier()')
    }

    # Check for missing values
    for (i in c(treatment, match.on)) {
      if (anyNA(dataset[[i]])) {
        customStop("missing values are not allowed in the treatment or matching variables.", 'makeFrontier()')
      }
    }

    # Check treatment
    if (length(unique(dataset[[treatment]])) != 2) {
        customStop('the treatment must be a binary variable (ideally 0/1).', 'makeFrontier()')
    }
}
