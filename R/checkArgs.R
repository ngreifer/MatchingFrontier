checkArgs <- function(QOI, metric, outcome = NULL, ...){
    if (!is.null(outcome)) {
        customWarning("the outcome argument is defunct and will be ignored. See ?estimateEffects.", 'makeFrontier()')
    }

    if (length(QOI) != 1 || !is.character(QOI)) {
        customStop("'QOI' must be a string.")
    }
    QOI <- try(match_arg(toupper(QOI), c("SATE", "FSATE", "SATT", "FSATT")), silent = TRUE)
    if (inherits(QOI, "try-error")){
        customStop("'QOI' must be either 'SATE', 'FSATE', 'SATT', or 'FSATT'.", 'makeFrontier()')
    }

    if (length(metric) != 1 || !is.character(metric)) {
        customStop("'metric' must be a string.")
    }
    metric <- try(match_arg(tolower(metric), c("l1", "l1median", "l2", "l2median", "mahal", "euclid", "custom", "energy")), silent = TRUE)
    if (inherits(metric, "try-error")){
        customStop("'metric' must be either 'l1', 'l1median', 'l2', 'l2median', 'mahal', 'euclid', 'custom', or 'energy'.", 'makeFrontier()')
    }

    acceptable_dist <- c("FSATE", "FSATT")
    acceptable_bin <- c("FSATE", "SATT")
    acceptable_energy <- c("SATE", "FSATE", "SATT")

    acceptable_combinations <- list(
        mahal = acceptable_dist,
        euclid = acceptable_dist,
        custom = acceptable_dist,
        l1 = acceptable_bin,
        l2 = acceptable_bin,
        l1median = acceptable_bin,
        l2median = acceptable_bin,
        energy = acceptable_energy
    )

    bad_combination <- !QOI %in% acceptable_combinations[[metric]]

    if (bad_combination) {
        msg <- paste0('the ', metric, ' theoretical frontier is only presently calculable for the ',
                      word_list(acceptable_combinations[[metric]]), '.')
        customStop(msg, 'makeFrontier()')
    }

}
