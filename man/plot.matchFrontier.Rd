\name{plot.matchFrontier}
\alias{plot.matchFrontier}

\title{
Plot the Balance-Sample Size Frontier
}
\description{
Plots the balance-sample size frontier for the supplied \code{matchFrontier} object. Can also plot covariate-specific mean differences and means to assess balance on individual covariates across the frontier.
}
\usage{
\method{plot}{matchFrontier}(x, covs = NULL, stat = "std-diff", n.estimated = 250, ...)
}

\arguments{
  \item{x}{
a \code{matchFrontier} object; the output of a call to \fun{makeFrontier}.
}
  \item{covs}{
the names of covariates for which mean differences or means are to be displayed. If \code{NULL} (the default), will display the frontier for the balance metric used in constructing the frontier.
}
  \item{stat}{
if \code{covs} is not \code{NULL}, the statistic to display for each variable across the frontier. Allowable options include \code{"std-diff"} for standardized mean differences (the default), \code{"diff"} for raw mean differences, \code{"ks"} for Kolmogorov-Smirnov statistics, \code{"std-mean"} for standardized means, \code{"mean"} for raw means, and \code{"ks-target"} for sample-target Kolmogorov-Smirnov statistics. See Details. Abbreviations allowed.
}
  \item{n.estimated}{
if \code{covs} is not \code{NULL}, the number of points along the frontier for which to compute balance statistics. Computing statistics at every point on the frontier can be time-consuming for large frontiers, especially when computing Kolmogorov-Smirnov statistics, so using values less than the number of points can save time at the expense of providing less detail in the frontier plot. The default is 250 or the number of points, whichever is smaller. Ignored when \code{covs = NULL}.
}
  \item{\dots}{
if \code{covs} is \code{NULL}, further arguments passed to \pkgfun2{ggplot2}{geom_line}{ggplot2::geom_line} to control the appearance of the plot (e.g., \code{color}, \code{size}, etc.). Ignored otherwise.
}
}
\details{
\code{plot.matchFrontier()} uses \pkgfun2{ggplot2}{geom_line}{ggplot2::geom_line} to display the balance-sample size frontier. The plot has two x-axes: the bottom one indicates the number of units dropped, and the top one indicates the number of units remaining. Which group this "number of units" refers to depends on the \code{QOI} and \code{metric} supplied to the original call to \fun{makeFrontier}.

When \code{covs} is \code{NULL}, the y-axis indicates the value of the imbalance metric at the given sample size; otherwise, it indicates the value of the balance statistic for the covariates supplied to \code{covs} depending on the arguments to \code{stat}.

When covariates are supplied to \code{covs}, one of six plots can be produced:
\enumerate{
\item standardized mean differences between treatment groups (\code{stat = "std-diff"}, the default).
\item raw mean differences between treatment groups (\code{stats = "diff"})
\item Kolmogorov-Smirnov statistics between the treatment groups (\code{stat = "ks"})
\item standardized means, which are equivalent to standardized mean differences between the full matched sample and the target sample (\code{stat = "std-mean"})
\item raw means in the matched sample (\code{stat = "mean"})
\item Sample-target Kolmogorov-Smirnov statistics, the Kolmogorov-Smirnov statistic between the full matched sample and the target sample (\code{stat = "ks-target"})
}

The target sample depends on the argument to \code{QOI} in the original call to \fun{makeFrontier}: when \code{"SATE"} or \code{"FSATE"}, the target sample is the full sample prior to matching; when \code{"SATT"} or \code{"FSATT"}, the target sample is the treated group prior to matching. The target sampel also determines the standardization factor when compute standardized means or mean differences.

When \code{stat = "mean"}, a dashed line will appear on the plot for each covariate indicating its mean in the target sample. Otherwise, an additional point will be plotted indicating the statistic computed prior to matching; for energy distance- and bin-based frontiers, this is the start of the frontier, but for pair distance-based frontier, the frontier starts after matching, so the additional point corresponds to balance prior to the initial point of the frontier.
}
\value{
A \code{ggplot} object that can be further manipulated using \pkg{ggplot2} functionality. For example, to zoom in on an area of the plot, \pkgfun{ggplot2}{coord_cartesian} can be used.
}

\seealso{
\fun{makeFrontier}, \fun{summary.matchFrontier}, \fun{plot.frontierEstimates}
}
\examples{
data("lalonde", package = "MatchIt")

#Pair distance frontier for FSATT
f1 <- makeFrontier(treat ~ age + educ + married + re74,
                   data = lalonde, QOI = "FSATT",
                   metric = "mahal")

#Plot the balance metric frontier
plot(f1)

#Plot standardized mean differences
plot(f1,
     covs = c("age", "educ", "married", "re74"),
     stat = "std-diff")

#Plot standardized deviations from the mean
#of the treated group in the original
#sample (prior to matching)
plot(f1,
     covs = c("age", "educ", "re74"),
     stat = "std-mean")

#Plot covariate means across the frontier;
#make sure variables are on the same scale!
plot(f1,
     covs = c("married", "nodegree"),
     stat = "mean")
}
