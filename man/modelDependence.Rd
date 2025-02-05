\name{modelDependence}
\alias{modelDependence}
\alias{print.modelDependenceBounds}
\alias{plot.modelDependenceBounds}

\title{
Compute Model Dependence Bounds
}
\description{
\code{modelDependence()} computes model dependence bounds for a single dataset using either the extreme bounds procedure or Athey-Imbens procedure. The width of the bounds indicates the degree of model dependence present in the dataset. Wide bounds imply that treatment effects could vary greatly depending on the outcome model used and therefore depend heavily on it.
}
\usage{
modelDependence(base.form = NULL,
                data,
                treatment,
                outcome,
                covariates,
                weights = NULL,
                method = c("extreme-bounds", "athey-imbens"),
                model.dependence.ests = 101,
                specifications = NULL,
                cutpoints = NULL,
                cutpoint.method = c("mean", "median", "segmented"),
                QOI = "FSATT",
                verbose = TRUE,
                seed = NULL)

\method{plot}{modelDependenceBounds}(x, ...)

\method{print}{modelDependenceBounds}(x, ...)
}
\arguments{
  \item{base.form}{
a \code{\link{formula}} relating the outcome to the treatment and covariates, e.g., \code{Y ~ A + X1 + X2 + ...} where \code{Y} represents the outcome variable, \code{A} the treatment variable, and \code{X1} and \code{X2} the covariates. This should correspond to the "default" or originally intended outcome model whose sensitivity is being assessed. If omitted, will just use a simple regression of the outcome on the treatment.
}
  \item{data}{
a data frame containing the variables named in \code{base.form}, \code{treatment}, \code{outcome}, and \code{covariates}.
}
  \item{treatment}{
the name of the treatment variable. If missing and \code{base.form} is specified and only has one right-hand side variable, that variable will be considered the treatment.
}
  \item{outcome}{
the name of the outcome variable. Ignored if \code{base.form} is supplied and required otherwise.
}
  \item{covariates}{
a character vector containing the names of the covariates over which model dependence is to be assessed. Can also be a one-sided formula with the covariates on the right-hand side. If omitted, all variables in \code{base.form} that are not the treatment or outcome will be considered covariates.
}
  \item{weights}{
an optional vector of weights (e.g., matching weights) to be included in the fitted models.
}
  \item{method}{
the method of computing model dependence bounds. Can be one of \code{"extreme-bounds"} or \code{"athey-imbens"} for the extreme bounds or Athey-Imbens procedures, respectively. Abbreviations allowed. Default is extreme bounds. See Details.
}
  \item{model.dependence.ests}{
when \code{method = "extreme-bounds"}, the number of different model specifications to use. Higher numbers give more accurate bounds but take longer. Exactly one of \code{model.dependence.ests} or \code{specifications} must be supplied.
}
  \item{specifications}{
when \code{method = "extreme-bounds"}, a list of \code{formula}s containing all the model specifications to consider in computing the extreme bounds. If omitted, these will be randomly generated using up to third degree polynomials and two-way interactions between covariates. Exactly one of \code{model.dependence.ests} or \code{specifications} must be supplied.
}
  \item{cutpoints}{
when \code{method = "athey-imbens"}, a named list of the cutpoints for each covariate that has more than two unique values. A cutpoint of \code{NA} means the unique values of the variable will be used. For any covariates without cutpoints specified, the cutpoints will computed using the method supplied to \code{cutpoint.method}.
}
  \item{cutpoint.method}{
when \code{method = "athey-imbens"}, the name of the method to create cutpoints for covariates without entries in \code{cutpoints}. Options include \code{"mean"} (the default) to use the covariate mean as the cutpoint, \code{"median"} to use the covariate median as the cutpoint, and \code{"segmented"} to estimate the cutpoints with segmented (piecewise) regression. Abbreviations allowed.
}
  \item{QOI}{
the quantity of interest (i.e., estimand); allowable options include \code{'SATE'}, \code{'FSATE'}, \code{'SATT'}, and \code{'FSATT'}. See \fun{makeFrontier}. Default is \code{'FSATT'}. This affects for whom treatment effects are estmated.
}
  \item{verbose}{
\code{logical}; whether to print additional information about the procedure.
}
  \item{seed}{
the seed for the random number generator. If \code{NULL}, the current \link[=set.seed]{seed} is used.
}
  \item{x}{
a \code{modelDependenceBounds} object, the output of a call to \code{modelDependence()}.
}
  \item{\dots}{
for \code{print.modelDependenceBounds}, further arguments passed to \fun{print}. Ignored otherwise.
}
}
\details{
The simplest way to estimate the model dependence bounds is to supply \code{base.form} with the treatment and covariates on the right hand side and \code{treatment} supplied.

The extreme bounds procedure (Leamer, 1983) works by fitting many different models of various functional forms with many combinations of the covariates. When \code{specifications} is not supplied, these models are randomly generated by appending to the model in \code{base.form} sets of squares, cubes, and two-way interactions of the covariates with \code{model.dependence.ests} controlling how many such models are estimated. More is always better but will take longer. The 5th and 95th percentiles of the treatment effects from all the specifications are used as the model dependence bounds.

The Athey-Imbens procedure (Athey & Imbens, 2015) works by splitting covariates at points determined by \code{cutpoints} or \code{cutpoint.method}, fitting the \code{base.form} model within each subset, averaging the two effect estimates, and computing the standard deviation of the resulting effect estimates. The base treatment effect estimate (fit in the full sample) minus and plus this standard deviation are used as the model dependence bounds. \code{cutpoints} can be supplied or estimated for each covariate. For covariates with exactly two unique values, the data will be split according to those values. For covariates with more than two values, if a cutpoint is supplied, it will be used to split the dataset; otherwise, a cutpoint will be estimated as the mean or median of the covariate or using segmented regression as implemented in \pkgfun2{segmented}{segmented}{segmented::segmented} depending on the argument to \code{cutpoint.method}.

It is important to consider the interpretation of \code{base.form}. If there is no confounding by the covariates (e.g., because of the result of a randomized experiment or successful matching), \code{base.form} can specified simply as \code{Y ~ A}, with model dependence referring to all the possible adjustment sets and their functional forms. If there is confounding and \code{base.form} is specified as \code{Y ~ A} (or with only a subset of the confounding covariates included in the model), then model dependence concerns different adjustment sets, and therefore muddles dependence on model specification and covariate choice, which is undesirable. \code{base.form} should therefore contain all the covariates required to eliminate confounding, with the model dependence referring to the possible additional terms that could be added to the model to change its functional form. Put simply, the adjustment set should always be in \code{base.form}, but after randomization or successful matching the adjustment set might be empty.
}
\value{
A \code{modelDependenceBounds} object, which is numeric vector with the lower and upper model dependence bounds and attributes concerning the fit. The attribute \code{"method"} will contain the method supplied to \code{method}. The attribute \code{base.est} will contain the effect estimated using \code{base.form} with no modification. The attribute \code{"estimates"} will contain the effect estimates used in computing the bounds. For the extreme bounds procedure, the attribute \code{"sepcifications"} contains a list of the model specifications used, each as a string representing the model formula. This can be supplied to future calls to \code{modelDependence()} or \fun{estimateEffects}. For the Athey-Imbens procedure, the attribute \code{"sigma"} contains the standard deviation of the effect estimates that, when added to and subtracted from the estimate resulting from the model specified in \code{base.form}, forms the bounds.

The \code{print()} method displays the bounds as well as their difference. Note that the difference is not actually stored in the \code{modelDependenceBounds} object but can be computed using \fun{diff}.

The \code{plot()} method displays a histogram of the effect estimates when \code{method = "extreme-bounds"} and a dot plot of the estimates when \code{method = "athey-imbens"}. Red vertical lines indicate the reported bounds, and the blue vertical line indicates the base estimate. The resulting object is a \code{ggplot} object that can be further manipulated using \pkg{ggplot2} functionality.
}
\references{
Athey, S., & Imbens, G. (2015). A measure of robustness to misspecification. \emph{American Economic Review}, 105(5), 476-80.

Leamer, E. E. (1983). Let's take the con out of econometrics. \emph{The American Economic Review}, 73(1), 31-43.
}
\seealso{
\fun{estimateEffects}, \CRANpkg{ExtremeBounds}
}
\examples{
library(MatchIt)
data("lalonde", package = "MatchIt")

#Model dependence prior to matching
b1 <- modelDependence(re78 ~ treat + age + educ + married +
                        race + nodegree + re74 + re75,
                      data = lalonde, treatment = "treat",
                      method = "extreme-bounds", seed = 12345)
b1 # diff = 823.8418
plot(b1)

#Model dependence after matching
m.out <- matchit(treat ~ age + educ + married + race +
                   nodegree + re74 + re75, data = lalonde,
                 distance = "mahalanobis", replace = TRUE)

b2 <- modelDependence(re78 ~ treat + age + educ + married +
                        race + nodegree + re74 + re75,
                      data = lalonde, treatment = "treat",
                      method = "extreme-bounds", seed = 12345,
                      weights = m.out$weights)
b2 # diff = 684.6681
plot(b2)

#Athey-Imbens bounds using various cutpoints
#  Here we use cutpoints of 100 for re74 and re75 and split
#  other covariates at their means
b3 <- modelDependence(re78 ~ treat + age + educ + married +
                        race + nodegree + re74 + re75,
                      data = lalonde, treatment = "treat",
                      method = "athey-imbens",
                      cutpoints = list(re74 = 100, re75 = 100),
                      cutpoint.method = "mean")
b3
plot(b3)
}
