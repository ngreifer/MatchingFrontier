\name{makeFrontier}
\alias{makeFrontier}
\alias{print.matchFrontier}
\alias{makeFrontier.formula}
\alias{makeFrontier.data.frame}

\title{
Compute the balance-sample size frontier
}
\description{
\code{makeFrontier()} computes the balance-sample size frontier, a function relating the size of a subset of the original sample to the least imbalance possible for that subset given an imbalance metric and a quantity of interest (i.e., estimand). The user can visualize the frontier using \code{plot()} methods, then pick a point along the frontier that represents what they consider to yield the optimal bias-variance tradeoff or estimate treatment effects at each point.
}
\usage{
makeFrontier(x, ...)

\method{makeFrontier}{formula}(formula,
             data,
             QOI = 'FSATT',
             metric = 'mahal',
             breaks = NULL,
             distance.mat = NULL,
             verbose = TRUE, ...)

\method{makeFrontier}{data.frame}(x,
             treatment,
             match.on,
             QOI = 'FSATT',
             metric = 'mahal',
             breaks = NULL,
             distance.mat = NULL,
             verbose = TRUE, ...)

\method{print}{matchFrontier}(x, ...)

}
\arguments{
  \item{formula}{
a two-sided \code{formula} object containing the treatment and covariates to be used in creating the distance measure used in the matching, e.g., \code{A ~ X1 + X2 + ...} where \code{A} represents the treatment variable and \code{X1} and \code{X2} are covariates.
}
  \item{x}{
for \code{makeFrontier()}, a data frame containing the variables named in \code{treatment} and \code{match.on}.

For \code{print()}, a \code{matchFrontier} object; the output of a call to \code{matchFrontier()}.
}
  \item{data}{
a data frame containing the variables included in \code{formula}.
}
  \item{treatment}{
a string containing the name of the treatment variable. Only required for the \code{data.frame} method.
}
  \item{match.on}{
a string containing the name of the covariates variable to match on. Only required for the \code{data.frame} method.
}
  \item{QOI}{
the quantity of interest (i.e., estimand); allowable options include \code{'SATE'}, \code{'FSATE'}, \code{'SATT'}, and \code{'FSATT'}. See Details. Default is \code{'FSATT'}. Not all \code{QOI}s can be used with all \code{metric}s.
}
  \item{metric}{
the imbalance metric to be minimized; allowable options include \code{'L1'}, \code{'L1median'}, \code{'L2'}, \code{'L2median'}, \code{'mahal'}, \code{'euclid'}, \code{'custom'}, and \code{'energy'}. See Details. Default is \code{'mahal'}. Not all \code{metric}s can be used with all \code{QOI}s.
}
  \item{breaks}{
when \code{metric} is \code{'L1'} or \code{'L2'}, a named list containing a single number giving the number of bins for each or a string corresponding to a method of computing the number of bins for each covariate; for example, \code{list(X1 = 4, X3 = "sturges")}. Allowable strings include \code{"sturges"}, \code{"scott"}, and \code{"fd"}, which use the functions \pkgfun2{grDevices}{nclass}{nclass.Sturges}, \pkgfun2{grDevices}{nclass}{nclass.scott}, and \pkgfun2{grDevices}{nclass}{nclass.FD}, respectively. The default is \code{"sturges"} for variables that are not listed or if no argument is supplied. Can also be a single value to be applied to all numeric variables. Similar to the \code{cutpoints} argument of \code{MatchIt::matchit()} with \code{\link[MatchIt:method_cem]{method = "cem"}}. Ignored for other \code{metric}s.
}
  \item{distance.mat}{
when \code{metric} is \code{'custom'} or \code{'energy'}, a distance matrix used to define the distance between units. For \code{metric = 'custom'}, should be a distance matrix with one row for each treated unit and one column for each control unit. For \code{metric = 'energy'}, should be a full N by N distance matrix with a row and column for each unit. Optional when \code{metric = 'energy'}; if not specified, the Euclidean distance computed on the scaled covariates will be used instead. See Details.
}
  \item{verbose}{
\code{logical}; whether to print additional information about the frontier constructing process.
}
  \item{\dots}{
arguments passed to other methods (ignored).
}
}
\details{
\code{makeFrontier()} computes the matching frontier for a given dataset, quantity of interest (QOI), and imbalance metric. The frontier describes the relationship between the size of the remaining sample and the lowest imbalance possible given that sample size. How the frontier is formed depends on the QOI and imbalance metric. Three types of frontier are currently available depending on the argument supplied to \code{metric}: pair distance-based frontiers (\code{'mahal'}, \code{'euclid'}, and \code{'custom'}), bin-based frontiers (\code{'L1'}, \code{'L1median'}, \code{'L2'}, and \code{'L2median'}) and energy distance-based frontiers (\code{'energy'}). These are described briefly below and in more detail \code{vignette("MatchingFrontier")}.

\subsection{Pair distance-based frontiers}{
With a pair distance-based frontier, each unit is paired with the closest unit of the opposite treatment group, where the distance is measured as the Mahalanobis distance, the Euclidean distance, or a user-supplied distance, as determined by the \code{metric} argument. Pairing is done with replacement so units can be paired to multiple units at a time. The average distance between each unit and its paired unit is used as the imbalance metric. The frontier is formed by dropping one unit at a time, starting with the unit that is farthest from its paired unit. If two units are equally far from their respective paired units, they are both dropped. In this way, the imbalance metric will decrease monotonically until no units remain or all remain pairs contain identical members.

These frontiers can only be used when the QOI is the FSATE or FSATT. The distance between units used to form the pairs and compute the imbalance metric can be computed the Mahalanobis distance (\code{metric = 'mahal'}), the Euclidean distance (\code{metric = 'euclid'}), or a user-supplied distance (\code{metric = 'custom'}). When (\code{metric = 'custom'}), a distance matrix must be supplied to \code{distance.mat} with a row for each treated unit and a column for each control unit.
}

\subsection{Bin-based frontiers}{
Bin-based frontiers involve assigning each unit into a bin based on its unique combination of coarsened covariates. The frontier is formed by dropping one unit at a time, in particular, the unit that when dropped would yield the greatest decrease in the imbalance metric. Bin-based frontiers are formed by subset selection rather than pairing; units are not paired with each other, and the initial point of the frontier is computed on the unmodified sample (in contrast to pair distance-based frontiers, where the initial point is computed after pairing). Bin-based frontiers can be thought of as a close relative to coarsened exact matching, as implemented in \pkg{cem} and \pkg{MatchIt} (i.e., using \code{matchit()} with \code{method = "cem"}).

These frontiers can be used when the QOI is the FSATE or SATT. The imbalance metrics include the L1 and L2 statistic, which summarize the discrepancy between the relative proportions of treated and control units in each bin. Bins can be computed by supplying an argument to \code{breaks}, or they can be computed automatically as the bins that yield the median value of the imbalance metric across a random sample of binnings. For the latter, \code{metric} should be specified as \code{'L1median'} or \code{'L2median'}. Because this involves a random process, a \link[=set.seed]{seed} should be set prior to running to ensure results are replicable.
}

\subsection{Energy distance-based frontiers}{
Like bin-based frontiers, energy distance-based frontiers involve subset selection rather than pairing. The frontier is formed by dropping one unit at a time, in particular, the unit that when dropped would yield the lowest value of the imbalance metric. The imbalance metric is the energy distance, a scalar measure of the difference between two multivariate cumulative density functions. In each sample, each unit's contribution to the energy distance is computed, and the unit with the greatest contribution is dropped. Unlike the other frontiers, the energy distance-based frontier is non-monotonic; there is often a point when discarding additional units increases rather than decreases the imbalance metric.

These frontiers can be used when the QOI is the SATE, FSATE, or SATT. For the SATE, the energy distance is computed as Huling and Mak's (2020) "improved" energy distance, which is the sum of the energy distance between the treated group and the original full sample, the energy distance between the control group and the full sample, and the energy distance between the treated and control groups. Units from either group may be dropped. As the energy distance decreases, the resulting groups each increasingly resemble the full sample and each other. For the FSATE, the energy distance is computed only between the treated and control groups. For the SATT, the energy distance as computed only between the treated and control groups, but only control units can be dropped.
}

\subsection{Quantities of interest (QOIs)}{
The \code{QOI} argument can be set to control the estimand, i.e., the group to which the effect is meant to generalize. The sample average treatment treatment effect (SATE) is the average effect of treatment in the entire sample. The feasible sample average treatment effect (FSATE) is the average effect of treatment in the sample remaining after units form both groups are dropped. The sample average treatment effect in the treated (SATT) is the average effect of treatment in the sample of those who were treated. The feasible sample average treatment effect in the treated (FSATT) is the average effect of treatment in the sample of treated units after some treated units are dropped.
}

}
\value{
A \code{matchFrontier} object, which contains the matching frontier as well as information about the the function call.

\item{frontier}{the matching frontier, which contains four components:
\itemize{
\item \code{drop.order}: a list containing the index of the units dropped at each point on the frontier. The first entry will be empty, indicating that no units have been dropped.
\item \code{Xs}: a vector containing the cumulative number of units dropped at each point on the frontier. The first value will be 0.
\item \code{Ys}: a vector containing the value of the imbalance metric at each point on the frontier.
\item \code{distance.mat}: when using a pair distance-based frontier or energy distance-based frontier, the distance matrix used to compute the imbalance statistic.
}}
\item{treatment}{the name of the treatment variable.}
\item{QOI}{the argument supplied to \code{QOI}.}
\item{metric}{the argument supplied to \code{metric}.}
\item{dataset}{the dataset containing the units.}
\item{match.on}{a vector of the names of the covariates used in the matching.}
\item{matched.to}{with a pair distance-based metric, the indices of the units to which each units is paired. When \code{QOI} is \code{'FSATT'}, control units will have \code{NA} values.}
\item{call}{the original call to \code{makeFrontier()}.}
\item{n}{the number of units present at the first point in the frontier (i.e., the original sample). For the FSATT, this is the number of treated units; for the SATT, this is the number of control units.}

}
\references{
King, G., Lucas, C., & Nielsen, R. A. (2017). The Balance-Sample Size Frontier in Matching Methods for Causal Inference. \emph{American Journal of Political Science}, 61(2), 473-489. \doi{10.1111/ajps.12272}

}

\seealso{
\fun{summary.matchFrontier} for viewing a summary of the frontier.

\fun{plot.matchFrontier} for plotting the frontier.

\fun{generateDataset} and \fun{frontier_to_matchit} for extracting a dataset or \code{matchit} object at one point on the frontier.

\fun{estimateEffects} for estimating effects along the frontier.

\pkgfun2{MatchIt}{matchit}{MatchIt::matchit} for a more traditional form of matching at a single point on the theoretical frontier.
}
\examples{
data("lalonde", package = "MatchIt")

#Pair distance frontier for FSATT
f1 <- makeFrontier(treat ~ age + educ + married + re74,
                   data = lalonde, QOI = "FSATT",
                   metric = "mahal")

f1
plot(f1)
summary(f1)
summary(f1, N = c(20, 40, 60))

#L1 distance frontier for FSATE
f2 <- makeFrontier(treat ~ age + educ + married + re74,
                   data = lalonde, QOI = "FSATE",
                   metric = "L1")

f2
plot(f2)
summary(f2)

#Energy distance frontier for SATE
f3 <- makeFrontier(treat ~ age + educ + married + re74,
                   data = lalonde, QOI = "SATE",
                   metric = "energy")

f3
plot(f3)
summary(f3)
summary(f3, Ndrop = seq(300, 500, by = 50))
}
