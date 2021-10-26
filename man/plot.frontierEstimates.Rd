\name{plot.frontierEstimates}
\alias{plot.frontierEstimates}

\title{
Plot effect estimates along the frontier
}
\description{
Plots the effect estimates along the frontier to see how the effect changes with different remaining sample sizes. The plot includes the effect estimates as well as confidence interval or model dependence bands. The bands are produced using \pkgfun2{ggplot2}{geom_ribbon}{ggplot2::geom_ribbon}.
}
\usage{
\method{plot}{frontierEstimates}(x,
     band, band.col = "red",
     band.border.col = "red",
     band.alpha = .5, ...)
}

\arguments{
  \item{x}{
a \code{frontierEstimates} object; the output of a call to \fun{estimateEffects}.
}
  \item{band}{
which type of band to display. Can be \code{"none"} for no band, \code{"confidence"} for confidence bands, or \code{"model-dependence"} for model dependence bands (only if requested in the original call to \fun{estimateEffects}). The default is \code{"model-dependence"} when model dependence bounds were requested in the original call to \code{estimateEffects()} and \code{"confidence"} otherwise.
}
  \item{band.col}{
the desired color of the band, if present. Passed to the \code{fill} argument of \pkgfun2{ggplot2}{geom_ribbon}{ggplot2::geom_ribbon}.
}
  \item{band.border.col}{
the desired color of the band borders if present. Passed to the \code{color} argument of \pkgfun2{ggplot2}{geom_ribbon}{ggplot2::geom_ribbon}.
}
  \item{band.alpha}{
the desired transparency of the band, if present. Passed to the \code{alpha} argument of \pkgfun2{ggplot2}{geom_ribbon}{ggplot2::geom_ribbon}.
}
  \item{\dots}{
further arguments passed to \pkgfun2{ggplot2}{geom_line}{ggplot2::geom_line} to control the appearance of the line (e.g., \code{color}, \code{size}, etc.).
}
}
\details{
\code{plot.frontierEstimates()} uses \pkgfun2{ggplot2}{geom_line}{ggplot2::geom_line} and \pkgfun2{ggplot2}{geom_ribbon}{ggplot2::geom_ribbon} to display the effect estimates along the frontier. The plot has two x-axes: the bottom one indicates the number of units dropped, and the top one indicates the number of units remaining. Which group this "number of units" refers to depends on the \code{QOI} and \code{metric} supplied to the original call to \fun{makeFrontier}. The y-axis indicates the effect estimate at the given sample size.

An additional point and error bar will be displayed corresponding to the effect and bounds estimated in the unadjusted sample (prior to matching or dropping units).
}
\value{
A \code{ggplot} object that can be further manipulated using \pkg{ggplot2} functionality.
}

\seealso{
\fun{plot.matchFrontier}, \fun{estimateEffects}
}
\examples{
#See examples at help("estimateEffects").
}
