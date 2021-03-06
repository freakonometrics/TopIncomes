% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/code_package_v4_juin_2019.R
\name{rpareto1}
\alias{rpareto1}
\title{Random generation of the Pareto 1 distribution}
\usage{
rpareto1(n, mu, alpha = 1/xi, xi = 1/alpha)
}
\arguments{
\item{n}{integer}

\item{mu}{numeric, lower bound, default is 1}

\item{alpha}{numeric, tail index}

\item{xi}{numeric, inverse of \code{alpha}}
}
\value{
generates \code{n} values of the Pareto 1 distribution
}
\description{
Random generation of the Pareto 1 distribution
}
\examples{
set.seed(123)
rpareto1(6, 1, alpha=1.5)
# expected value of a Pareto1(1,1.5)
mean(rpareto1(1e6, 1, alpha=1.5))
# theoretical expected value
1*(1.5)/(1.5-1)
}
\seealso{
\code{\link{dpareto1}}, \code{\link{ppareto1}} and \code{\link{qpareto1}}
}
