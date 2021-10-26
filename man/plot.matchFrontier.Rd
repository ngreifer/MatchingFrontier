\name{plot.matchFrontier}
\alias{plot.matchFrontier}

\title{
Plot the Balance-Sample Size Frontier
}
\description{
Plots the balance-sample size frontier for the supplied \code{matchFrontier} object. Can also plot covariate-specific mean differences and means to assess balance on individual covariates across the frontier.
}
\usage{
\method{plot}{matchFrontier}(x, covs = NULL, diffs = TRUE,
                   std = TRUE, ...)
}

\arguments{
  \item{x}{
a \code{matchFrontier} object; the output of a call to \fun{makeFrontier}.
}
  \item{covs}{
the names of covariates for which mean differences or means are to be displayed. If \code{NULL} (the default), will display the frontier for the balance metric used in constructing the frontier.
}
  \item{diffs}{
if \code{covs} is not \code{NULL}, whether to display mean differences for the covariates (\code{TRUE}) or the covariate means (\code{FALSE}).
}
  \item{std}{
if \code{covs} is not \code{NULL}, whether to display standardized statistics (i.e., standardized mean differences or standardized means) (\code{TRUE}) or raw statistics (i.e., mean differences or means) (\code{FALSE}). Note that it rarely makes sense to plot multiple variables on different scales without standardizing them.
}
  \item{\dots}{
if \code{covs} is \code{NULL}, further arguments passed to \pkgfun2{ggplot2}{geom_line}{ggplot2::geom_line} to control the appearance of the plot (e.g., \code{color}, \code{size}, etc.). Ignored otherwise.
}
}
\details{
\code{plot.matchFrontier()} uses \pkgfun2{ggplot2}{geom_line}{ggplot2::geom_line} to display the balance-sample size frontier. The plot has two x-axes: the bottom one indicates the number of units dropped, and the top one indicates the number of units remaining. Which group this "number of units" refers to depends on the \code{QOI} and \code{metric} supplied to the original call to \fun{makeFrontier}.

When \code{covs} is \code{NULL}, the y-axis indicates the value of the imbalance metric at the given sample size; otherwise, it indicates the (standardized) means or mean differences of the covariates supplied to \code{covs} depending on the arguments to \code{diffs} and \code{std}.

When covariates are supplied to \code{covs}, one of four plots can be produced:
\enumerate{
\item standardized mean differences between treatment groups (\code{diffs = TRUE, std = TRUE}, the default)
\item raw mean differences between treatment groups (\code{diffs = TRUE, std = FALSE})
\item standardized mean differences between the matched sample and the target sample (\code{diffs = FALSE, std = TRUE})
\item raw means in the matched sample (\code{diffs = FALSE, std = FALSE})
}

The target sample depends on the argument to \code{QOI} in the original call to \fun{makeFrontier}: when \code{"SATE"} or \code{"FSATE"}, the target sample is the full sample prior to matching; when \code{"SATT"} or \code{"FSATT"}, the target sample is the treated group prior to matching. The standardization factor used when \code{std = TRUE} is the standard deviation of the covariate in the target sample.

When displaying (standardized) mean differences between treatment groups, an additional point will be plotted indicating the statistic computed prior to matching; for energy distance- and bin-based frontiers, this is the start of the frontier, but for pair distance-based frontier, the frontier starts after matching, so the additional point corresponds to balance prior to the initial point of the frontier.

When displaying covariate means, a dashed line will appear on the plot for each covariate indicating its mean in the target sample. When \code{std = TRUE}, the standardized means are all zero in the target population by definition; deviations from zero indicate deviations from the target population mean in standard deviation units. As with mean differences, the target sample covariate means may not coincide with the start of the frontier for pair distance-based frontiers.
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
                   metric = "Mahal")

#Plot the balance metric frontier
plot(f1)

#Plot standardized mean differences
plot(f1,
     covs = c("age", "educ", "married", "re74"),
     diffs = TRUE, std = TRUE)

#Plot standardized deviations from the mean
#of the treated group in the original
#sample (prior to matching)
plot(f1,
     covs = c("age", "educ", "re74"),
     diffs = FALSE, std = TRUE)

#Plot covariate means across the frontier;
#make sure variables are on the same scale!
plot(f1,
     covs = c("married", "nodegree"),
     diffs = FALSE, std = FALSE)
}
