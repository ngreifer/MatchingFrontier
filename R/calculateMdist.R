calculateMdist <- function(covs.mat, treat.vec){
    ## calculate the inverse covariance matrix

    icv <- generalized_inverse(cov(covs.mat))

    treated.ind <- which(treat.vec == 1)
    control.ind <- which(treat.vec == 0)

    Ts <- covs.mat[treated.ind,, drop = FALSE]
    Cs <- covs.mat[control.ind,, drop = FALSE]

    d <- matrix(NA_real_, nrow = length(treated.ind), ncol = length(control.ind),
                dimnames = list(rownames(covs.mat)[treated.ind],
                                rownames(covs.mat)[control.ind]))

    for (t in 1:nrow(Ts)) {
        d[t,] <- sqrt(mahalanobis(Cs, Ts[t,], icv, inverted = TRUE))
    }

    return(d)
}
