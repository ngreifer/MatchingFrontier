checkArgs <- function(QOI, metric, outcome = NULL, ...){
    if (!is.null(outcome)) {
        customWarning("the outcome argument is defunct and will be ignored. See ?estimateEffects.", 'makeFrontier()')
    }

    QOI <- try(match_arg(QOI, c("SATE", "FSATE", "SATT", "FSATT")), silent = TRUE)
    if (inherits(QOI, "try-error")){
        customStop("QOI must be either 'SATE', 'FSATE', 'SATT', or 'FSATT'.", 'makeFrontier()')
    }

    metric <- try(match_arg(metric, c("L1", "L1median", "L2", "L2median", "Mahal", "Euclid", "Custom", "Energy")), silent = TRUE)
    if (inherits(metric, "try-error")){
        customStop("'metric' must be either 'L1', 'L1median', 'L2', 'L2median', 'Mahal', 'Euclid', 'Custom', or 'Energy'.", 'makeFrontier()')
    }

    acceptable_dist <- c("FSATE", "FSATT")
    acceptable_bin <- c("FSATE", "SATT")
    acceptable_energy <- c("SATE", "FSATE", "SATT")
    acceptable_combinations <- list(
        Mahal = acceptable_dist,
        Euclid = acceptable_dist,
        Custom = acceptable_dist,
        L1 = acceptable_bin,
        L2 = acceptable_bin,
        L1median = acceptable_bin,
        L2median = acceptable_bin,
        Energy = acceptable_energy
    )

    bad_combination <- !QOI %in% acceptable_combinations[[metric]]

    if (bad_combination) {
        msg <- paste0('the ', metric, ' theoretical frontier is only presently calculable for the ',
                      word_list(acceptable_combinations[[metric]]), '.')
        customStop(msg, 'makeFrontier()')
    }

}
