\name{summary.matchFrontier}
\alias{summary.matchFrontier}
\alias{print.summary.matchFrontier}

\title{
Summarize a \code{matchFrontier} object
}
\description{
\code{summary()} displays the sample sizes and imbalance metrics at the beginning and end of the frontier as well as at the point at which the imbalance metric is lowest.
}
\usage{
\method{summary}{matchFrontier}(object, ...)

\method{print}{summary.matchFrontier}(x, ...)
}

\arguments{
  \item{object}{
a \code{matchFrontier} object; the output of a call to \fun{makeFrontier}.
}
  \item{x}{
a \code{summary.matchFrontier} object; the output of a call to \code{summary.matchFrontier()}.
}
  \item{\dots}{
ignored.
}
}
\details{
The lowest value of the imbalance metric may not occur at the end of the frontier (i.e., with the smallest remaining sample size). \code{summary()} can be used to complement \fun{plot.matchFrontier} to identify which sample size yields the least imbalance.

A carrot (i.e., \code{^}) is printed below the column the corresponds to the \code{n} component of the \code{matchFrontier} object. It is this value that the \code{N} and \code{Ndrop} arguments refer to in \fun{frontier_to_matchit} and \fun{generateDataset}. For example, if the carrot is pointing at the \code{N control} column, then supplying a value of 50 to \code{N} in \fun{generateDataset} will produce a dataset with 50 control units (and however many treated units are present).
}
\value{
A \code{summary.matchFrontier} object, which contains the following components:

\item{Ntreated}{the number of treated units at the start, end, and best point on the frontier (best meaning with the lowest imbalance).}
\item{Ncontrol}{the number of control units at the start, end, and best point on the frontier.}
\item{N}{the total number of units at the start, end, and best point on the frontier.}
\item{Stat}{the value of the imbalance metric at the start, end, and best point on the frontier.}
\item{bestind}{the index of the best point on the frontier.}
\item{QOI}{the quantity of interest; the argument to \code{QOI} supplied in the original call to \fun{makeFrontier}.}
\item{metric}{the imbalance metric; the argument to \code{metric} supplied in the original call to \fun{makeFrontier}.}

}

\seealso{
\fun{makeFrontier}, \fun{plot.matchFrontier}
}
\examples{
#See examples at help("makeFrontier")
}
