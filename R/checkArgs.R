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

    acceptable_dist <- c("FSATE", "FSATT")
    acceptable_bin <- c("FSATE", "SATT")
    acceptable_energy <- c("SATE", "FSATE", "SATT")

    acceptable_combinations <- list(
        dist = acceptable_dist,
        l1 = acceptable_bin,
        l2 = acceptable_bin,
        energy = acceptable_energy
    )

    bad_combination <- !QOI %in% acceptable_combinations[[metric]]

    if (bad_combination) {
        msg <- sprintf('%s frontiers are only presently calculable for the %s.',
                      switch(metric, "dist" = "pair distance-based",
                             "energy" = "energy distance-based", "bin-based"),
                             word_list(acceptable_combinations[[metric]]))
        customStop(msg, 'makeFrontier()')
    }

}
