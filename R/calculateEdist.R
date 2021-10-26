calculateEdist <- function(covs.mat, treat.vec){

    treated.ind <- which(treat.vec == 1)
    control.ind <- which(treat.vec == 0)

    Ts <- covs.mat[treated.ind,, drop = FALSE]
    Cs_t <- t(covs.mat[control.ind,, drop = FALSE])

    d <- matrix(NA_real_, nrow = length(treated.ind), ncol = length(control.ind),
                dimnames = list(rownames(covs.mat)[treated.ind],
                                rownames(covs.mat)[control.ind]))
    for (t in 1:nrow(Ts)) {
        d[t,] <- sqrt(colSums((Ts[t,] - Cs_t)^2))
    }
    return(d)
}
